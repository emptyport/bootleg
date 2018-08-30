'use strict';

let mgf = require('js-mgf');
let fastaParser = require('fasta-js');
let peptideCutter = require('peptide-cutter');
let pepMod = require('peptide-modifier');
let pepFrag = require('peptide-fragmenter');

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
  for(var i=0; i<sequenceList.length; i++) {
    let residue = sequenceList[i];
    mass += residue_masses[residue];
  }

  for (var i=0; i<modifications.length; i++) {
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

const parseSpectra = (data) => {
  let spectra = mgf.parse(data);
  for(var i=0; i<spectra.length; i++) {
    let charge = parseFloat(spectra[i].charge);
    let pepmass = parseFloat(spectra[i].pepmass);
    let mass = charge * pepmass - charge;
    spectra[i].neutral_mass = mass;
    let intensitySum = spectra[i].intensity.reduce((total, num) => {
      return total + num;
    });
    for(var j=0; j<spectra[i].intensity.length; j++) {
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
    let matchTemplate = {
      bestScore: 0,
      fdr: 1,
      bestMatch: 0,
      matches: []
    }
    this.matches.fill(matchTemplate);
  }

  addMatch(index, modType, missedCleavages, fastaIndex, score, mzError) {
    let newMatch = {
      modType: modType,
      missedCleavages: missedCleavages,
      fastaIndex: fastaIndex,
      score: score,
      mzError: mzError
    };
    let spectrumMatch = this.matches[index];
    if(score > spectrumMatch.bestScore) {
      this.matches[index].bestScore = score;
      this.matches[index].bestMatch = this.matches[index].matches.length;
    }
    this.matches[index].matches.push(newMatch);
  }

}

class Database {
  constructor(config) {
    this.minimallyModified = {};
    this.modified = {};
    this.maxMissedCleavages = config.missedCleavages;
    let missed = 0;
    while (missed <= config.missedCleavages) {
      this.minimallyModified[missed] = [];
      this.modified[missed] = [];
      missed++;
    }
  }

  addPeptide(peptide) {
    if(peptide.modifications.length === 0) {
      this.minimallyModified[peptide.missedCleavages].push(peptide);
    }
    else {
      let isMinimallyModified = true;
      peptide.modifications.map((mod) => {
        if(mod.type === "variable" && mod.name.toLowerCase() !== "oxidation") { isMinimallyModified = false; }
      });
      if(isMinimallyModified) {
        this.minimallyModified[peptide.missedCleavages].push(peptide);
      }
      else {
        this.modified[peptide.missedCleavages].push(peptide);
      }
    }
  }

  sortPeptides() {
    let missed = 0;
    while (missed <= this.maxMissedCleavages) {
      this.minimallyModified[missed] = sortByKey(this.minimallyModified[missed], 'mass');
      this.modified[missed] = sortByKey(this.modified[missed], 'mass');
      missed++;
    }
  }
}

const processPeptides = (accession, sequence, config, database) => {
  let cuttingOptions = {
    'enzyme': config.enzyme,
    'num_missed_cleavages': config.missedCleavages,
    'min_length': config.minLength,
    'max_length': config.maxLength
  }
  let cutter = new peptideCutter(cuttingOptions);

  let cleavedPeptides = cutter.cleave(sequence);
  for(var i=0; i<cleavedPeptides.length; i++) {
    var pep = cleavedPeptides[i];
    let modifications = pepMod.modify(pep.sequence, config.modifications, config.numVariableMods);
    for(var j=0; j<modifications.length; j++) {
      let mass = calculateMass(pep.sequence, modifications[j]);
      if(mass === -1) { continue; }
      if(mass < config.minMass || mass > config.maxMass) { continue; }
      let dbEntry = {
        accession: accession,
        sequence: pep.sequence,
        start: pep.start,
        end: pep.end,
        missedCleavages: pep.missed,
        modifications: modifications[j],
        mass: mass,
        fragments: pepFrag.fragment(pep.sequence, config.ionTypes, config.fragmentCharge, modifications[j])
      };
      database.addPeptide(dbEntry);
    }
  }
}

const parseFasta = (data, config) => {
  let fasta = new fastaParser(config.fastaOptions);
  let sequenceDB = fasta.parse(data);
  let l = sequenceDB.length;

  // This database is going to be split into minimally
  // modified peptides and modified peptides as well as
  // by the number of missed cleavages. The idea is to 
  // build in a cascading search to help with search speed
  console.log('Creating database...');
  let count = 0;
  let numIntervalsTotal = 20;
  let numIntervalsElapsed = 0;
  let interval = parseInt(l/numIntervalsTotal);
  let database = new Database(config);
  sequenceDB.map((dbEntry) => {
    if(count%interval===0) {
      let progress = numIntervalsElapsed * (100/numIntervalsTotal);
      process.stdout.write(progress+"% ");
      numIntervalsElapsed++;
    }
    count++;
    let entryAccession = dbEntry.accession;
    let entrySequence = dbEntry.sequence;
    let peptidesToAdd = processPeptides(entryAccession, entrySequence, config, database);
    if(config.generateDecoy) {
      let decoy_accession = dbEntry.accession + config.decoyTag;
      let decoy_sequence = dbEntry.sequence.split("").reverse().join("");
      let decoyPeptidesToAdd = processPeptides(decoy_accession, decoy_sequence, config, database);
    }
  });
  console.log('\n');
  sequenceDB = [];
  console.log('Sorting database...');
  database.sortPeptides();
  return database;  

}

const matchFragments = (mz, fragments, config)=> {
  let matches = [];
  for(var i=0; i<fragments.length; i++) {
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

const runCascade = (fastaDB, spectra, config, matches) => {
  let peptideList;
  for(let i=0; i<=config.missedCleavages; i++) {
    console.log('Running search against minimally modified peptides with '+i+' missed cleavages');
    peptideList = fastaDB.minimallyModified[i];
    runSearch(peptideList, spectra, config, matches);
    console.log('Running search against modified peptides with '+i+' missed cleavages');
    peptideList = fastaDB.modified[i];
    runSearch(peptideList, spectra, config, matches);
  }
}

const runSearch = (peptideList, spectra, config, matches) => {
  if(peptideList.length === 0) {
    console.log("N/A");
    return;
  }

  let startIdx = 0;
  let endIdx = 0;
  for(let i=0; i<spectra.length; i++) {
    let spectrum = spectra[i];
    let errorRange = config.precursorTol / 1000000 * spectrum.neutral_mass;
    let lowerBound = spectrum.neutral_mass - errorRange;
    let upperBound = spectrum.neutral_mass + errorRange;

    while(peptideList[startIdx].mass < lowerBound && startIdx<peptideList.length) {
      startIdx++;
    }
    while(peptideList[endIdx].mass <= upperBound && endIdx<peptideList.length) {
      endIdx++;
    }

    for(let j=startIdx; j<=endIdx; j++) {
      console.log(peptideList[j]);
    }
  }


}

module.exports.search = (spectra_data, fasta_data, config) => {

  // Let's start off by parsing the spectra and then sorting everything
  let spectra = parseSpectra(spectra_data);
  let spectraLength = spectra.length;
  let spectraMassArray = new Array(spectraLength);
  for(let i=0; i<spectraLength; i++) {
    spectraMassArray[i] = (spectra[i].neutral_mass);
  }
  config.minMass = spectraMassArray[0].neutral_mass;
  config.maxMass = spectraMassArray[spectraLength - 1].neutral_mass;

  // Now let's parse the fasta file
  let fastaDB = parseFasta(fasta_data, config);
  
  console.log('Matching peptides to spectra...');
  let matches = new Matches(spectraLength);
  runCascade(fastaDB, spectra, config, matches);
  



}

