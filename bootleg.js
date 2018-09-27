'use strict';

let mgf = require('js-mgf');
let pepFrag = require('peptide-fragmenter');
let mysql = require('sync-mysql');

const calculateMass = (sequence, modifications) => {
  const residue_masses = {
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "D": 115.02694,
    "C": 103.00919,
    "E": 129.04259,
    "Q": 128.05858,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841
  };

  const ILLEGAL = [ 'B', 'J', 'O', 'U', 'X', 'Z' ];

  let sequenceList = sequence.split("");

  // For now we don't handle ambiguous symbols, just return nothing
  let overlap = sequenceList.filter(val => -1 !== ILLEGAL.indexOf(val));
  if(overlap.length > 0) {
    return -1;
  }

  let mass = 18.010565;
  for(var i=0, sequenceListLength=sequenceList.length; i<sequenceListLength; i++) {
    let residue = sequenceList[i];
    mass += residue_masses[residue];
  }

  for (var i=0, modificationsLength=modifications.length; i<modificationsLength; i++) {
    let modMass = modifications[i].mass;
    mass += modMass;
  }

  return mass;
  
}

// https://stackoverflow.com/questions/8837454/sort-array-of-objects-by-single-key-with-date-value
const sortByKey = (array, key) => {
  return array.sort(function(a, b) {
      var x = a[key]; var y = b[key];
      return ((x < y) ? -1 : ((x > y) ? 1 : 0));
  });
}

const sortByKeyReverse = (array, key) => {
  return array.sort(function(a, b) {
      var x = a[key]; var y = b[key];
      return ((x > y) ? -1 : ((x < y) ? 1 : 0));
  });
}

const parseSpectra = (data) => {
  let spectra = mgf.parse(data);
  for(var i=0, spectraLength=spectra.length; i<spectraLength; i++) {
    let charge = parseFloat(spectra[i].charge);
    let pepmass = parseFloat(spectra[i].pepmass);
    let mass = charge * pepmass - charge;
    spectra[i].neutral_mass = mass;
    let intensitySum = spectra[i].intensity.reduce((total, num) => {
      return total + num;
    });
    for(var j=0, spectraIntensityLength=spectra[i].intensity.length; j<spectraIntensityLength; j++) {
      spectra[i].intensity[j] = spectra[i].intensity[j] / intensitySum;
    }
  }
  console.log('Sorting spectra...');
  let sortedSpectra = sortByKey(spectra, 'neutral_mass');
  return sortedSpectra;
}

class Matches {
  constructor(length) {
    this.matches = new Array(length);
    // Tried using the Array.fill method but it passes in the object
    // reference instead of a new object. The below method keeps each
    // array entry as it's own thing.
    for(let i=0; i<length; i++) {
      this.matches[i] = {
        bestScore: 0,
        fdr: 1,
        bestMatch: 0,
        bestIsDecoy: false,
        matchedPSMs: []
      };
    }
  }

  addMatch(index, currentPeptide, score, isDecoy) {
    let newMatch = {
      missedCleavages: currentPeptide.missedCleavages,
      score: score,
      peptide: currentPeptide
    };
    let spectrumMatch = this.matches[index];
    if(score > spectrumMatch.bestScore) {
      this.matches[index].bestScore = score;
      this.matches[index].bestMatch = this.matches[index].matchedPSMs.length;
      this.matches[index].bestIsDecoy = isDecoy;
    }
    this.matches[index].matchedPSMs.push(newMatch);
  }

}

const prepareMatches = (spectra, matches) => {
  for(let i=0, matchesLength=matches.matches.length; i<matchesLength; i++) {
    let match = matches.matches[i];
    if(match.matchedPSMs.length === 0) {
      matches.matches[i].spectrum = '';
      matches.matches[i].accession = '';
      matches.matches[i].sequence = '';
      matches.matches[i].modifications = '';
      matches.matches[i].bestScore = '';
      continue;
    }
    let modType = match.matchedPSMs[match.bestMatch].modType;
    let missedCleavages = match.matchedPSMs[match.bestMatch].missedCleavages;
    let peptide = match.matchedPSMs[match.bestMatch].peptide;
    matches.matches[i].spectrum = spectra[i].title;
    matches.matches[i].accession = peptide.accession;
    matches.matches[i].sequence = peptide.sequence;
    matches.matches[i].modifications = peptide.modifications;
  }
}


const matchFragments = (mz, fragments, config)=> {
  let matches = [];
  for(var i=0, fragmentsLength=fragments.length; i<fragmentsLength; i++) {
    let closestMatch = closestIdx(fragments[i], mz);
    if(Math.abs(mz[closestMatch] - fragments[i]) < config.fragmentTol) {
      matches.push(closestMatch);
    }
  }
  return matches;
}

// https://stackoverflow.com/questions/8584902/get-closest-number-out-of-array
const closestIdx = (num, arr) => {
  var mid;
  var lo = 0;
  var hi = arr.length - 1;
  while (hi - lo > 1) {
      mid = Math.floor ((lo + hi) / 2);
      if (arr[mid] < num) {
          lo = mid;
      } else {
          hi = mid;
      }
  }
  if (num - arr[lo] <= arr[hi] - num) {
      return lo;
  }
  return hi;
}

const calculateFDR = (matches, config) => {
  let l = matches.matches.length;
  let calculationArr = new Array(l);
  for(let i=0; i<l; i++) {
    calculationArr[i] = {
      score: matches.matches[i].bestScore,
      fdr: matches.matches[i].fdr,
      isDecoy: matches.matches[i].bestIsDecoy,
      idx: i
    }
  }
  let sortedArr = sortByKey(calculationArr, 'score');
  let fdrMask = sortedArr.map((entry) => {
    return entry.isDecoy;
  });

  let numDecoys = fdrMask.reduce((total, num) => {
    return total + num;
  }, 0);
  let currentFDR = numDecoys / l;

  let idxCutoff = 0;
  while(currentFDR > config.fdrCutoff) {
    numDecoys -= fdrMask[idxCutoff];
    l--;
    currentFDR = numDecoys / l;
    idxCutoff++;
  }

  for(let i=idxCutoff, sortedArrLength=sortedArr.length; i<sortedArrLength; i++) {
    let idxToSet = sortedArr[i].idx;
    matches.matches[idxToSet].fdr = 0;
  }
  console.log('Retaining '+l+' matches');
}

// https://medium.freecodecamp.org/how-to-factorialize-a-number-in-javascript-9263c89a4b38
const factorialize = (num) => {
  if (num < 0) 
        return -1;
  else if (num == 0) 
      return 1;
  else {
      return (num * factorialize(num - 1));
  }
}

const runSearch = (spectra, config, matches, connection) => {
  console.log(`Searching ${spectra.length} spectra...`);
  let count = 0;
  let numIntervalsTotal = 20;
  let numIntervalsElapsed = 0;
  let interval = parseInt(spectra.length/numIntervalsTotal);

  for(let i=0, spectraLength=spectra.length; i<spectraLength; i++) {
    if(count%interval===0) {
      let progress = numIntervalsElapsed * (100/numIntervalsTotal);
      process.stdout.write(progress+"% ");
      numIntervalsElapsed++;
    }
    count++;

    let spectrum = spectra[i];
    let errorRange = config.precursorTol / 1000000 * spectrum.neutral_mass;
    let lowerBound = spectrum.neutral_mass - errorRange;
    let upperBound = spectrum.neutral_mass + errorRange;

    let query = 
      `SELECT * FROM ${config.table_name}
      WHERE mass>=${lowerBound} AND mass<=${upperBound}
      AND missed_cleavages<=${config.missedCleavages}`;

    let peptideList = connection.query(query);

    for(let j=0; j<peptideList.length; j++) {
      let currentPeptide = peptideList[j];
      currentPeptide.fragments = pepFrag.fragment(currentPeptide.sequence, ['b','y'], [1], currentPeptide.modifications);

      let matchedFragments = {};
      // This is ugly looking. Iterate over ion types and then charge states
      for(var ionType in currentPeptide.fragments) {
        if(currentPeptide.fragments.hasOwnProperty(ionType)) {

          matchedFragments[ionType] = {};
          for(var chargeState in currentPeptide.fragments[ionType]) {
            if(currentPeptide.fragments[ionType].hasOwnProperty(chargeState)) {

              let fragmentMatches = matchFragments(spectrum.mz, currentPeptide.fragments[ionType][chargeState], config);
              matchedFragments[ionType][chargeState] = fragmentMatches;
            }
          }
        }
      }
      // This is essentially Morpheus scoring. We iterate over the matched
      // indices and count how many we have and add it to the score. Then
      // we add in the sum of the intensities (which are already normalized)
      // so we get the fraction of matched peak intensity
      let score = 0;
      for(var ionType in matchedFragments) {
        if(matchedFragments.hasOwnProperty(ionType)) {

          for(var chargeState in matchedFragments[ionType]) {
            if(matchedFragments[ionType].hasOwnProperty(chargeState)) {

              let matchIdxArr = matchedFragments[ionType][chargeState];
              let numIons = matchIdxArr.length;
              score += numIons;
              for(let k=0, matchIdxArrLength=matchIdxArr.length; k<matchIdxArrLength; k++) {
                score += spectrum.intensity[k];
              }
            }
          }
        }
      }
      let isDecoy = false;
      if(currentPeptide.accession.includes(config.decoyTag)) {
        isDecoy = true;
      }
      if(score > 0) {
        matches.addMatch(i, currentPeptide, score, isDecoy);
      }
    }
  }
  process.stdout.write("\n");

}

module.exports.search = (spectra_data, config) => {

  // Let's start off by parsing the spectra and then sorting everything
  let spectra = parseSpectra(spectra_data);
  let spectraLength = spectra.length;
  let spectraMassArray = new Array(spectraLength);
  for(let i=0; i<spectraLength; i++) {
    spectraMassArray[i] = (spectra[i].neutral_mass);
  }
  config.minMass = spectraMassArray[0].neutral_mass;
  config.maxMass = spectraMassArray[spectraLength - 1].neutral_mass;

  console.log('Matching peptides to spectra...');
  let matches = new Matches(spectraLength);

  let connection = new mysql({
    host     : config.server_address,
    port     : config.server_port,
    user     : config.server_username,
    password : config.server_password,
    database : config.db_name
  });

  runSearch(spectra, config, matches, connection);
  calculateFDR(matches, config);
  prepareMatches(spectra, matches);
  return(matches);
}

