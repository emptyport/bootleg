let config = require('./config.json');
let fs = require('fs');
let fastaParser = require('fasta-js');
let peptideCutter = require('peptide-cutter');
let pepMod = require('peptide-modifier');
let pepFrag = require('peptide-fragmenter');
let mgf = require('js-mgf');
const createCsvWriter = require('csv-writer').createObjectCsvWriter;

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

let fastaOptions = {
  'definition': 'gi|accession|description',
  'delimiter': '|'
};

let cuttingOptions = {
  'enzyme': 'trypsin',
  'num_missed_cleavages': config.missedCleavages,
  'min_length': config.minLength,
  'max_length': config.maxLength
}

let fasta = new fastaParser(fastaOptions);
let cutter = new peptideCutter(cuttingOptions);

fs.readFile('chicken.mgf', 'utf-8', function(err, data) {
  if (err) throw err;
  console.log('Read in mgf file');
  console.log('Analyzing mgf file...');
  let spectra = mgf.parse(data);
  let minMass = 10000;
  let maxMass = 0;
  for(var i=0; i<spectra.length; i++) {
    let charge = parseFloat(spectra[i].charge);
    let pepmass = parseFloat(spectra[i].pepmass);
    let mass = charge * pepmass - charge;
    spectra[i].neutral_mass = mass;
    spectra[i].matches = [];
    if(mass < minMass) { minMass = mass; }
    if(mass > maxMass) { maxMass = mass; }
    let intensitySum = spectra[i].intensity.reduce(getSum);
    for(var j=0; j<spectra[i].intensity.length; j++) {
      spectra[i].intensity[j] = spectra[i].intensity[j] / intensitySum;
    }
  }

  config.minMass = minMass;
  config.maxMass = maxMass;

  sortedSpectra = sortByKey(spectra, 'neutral_mass');
  let spectraLength = sortedSpectra.length;
  let spectraMassArray = new Array(spectraLength);
  for(var i=0; i<spectraLength; i++) {
    spectraMassArray[i] = (sortedSpectra[i].neutral_mass);
  }

  fs.readFile('chicken.fasta', 'utf-8', function(err, data){
    if(err) throw err;
    console.log('Read in fasta file');
    console.log('Searching...');

    let sequenceDB = fasta.parse(data);
    let l = sequenceDB.length;
    let interval = parseInt(l/50);
    for(var i=0; i<l; i++) {
      if(i%interval === 0) { process.stdout.write("#"); }
      let dbEntry = sequenceDB[i];
      let entry_accession = dbEntry.accession;
      let entry_sequence = dbEntry.sequence;
      let decoy_accession = dbEntry.accession + "_REVERSED";
      let decoy_sequence = dbEntry.sequence.split("").reverse().join("");
      
      let peptidesToSearch = [];
      peptidesToSearch.extend(processPeptides(entry_accession, entry_sequence, config));
      peptidesToSearch.extend(processPeptides(decoy_accession, decoy_sequence, config));

      for(var j=0; j<peptidesToSearch.length; j++) {
        let p = peptidesToSearch[j];
        if(p.mass<spectraMassArray[0] || p.mass>spectraMassArray[spectraLength]) { continue; }
        let errorRange = config.precursorTol / 1000000 * p.mass;
        let startMass = p.mass - errorRange;
        let endMass = p.mass + errorRange;
        let startIdx = closestIdx(startMass, spectraMassArray);
        let endIdx = closestIdx(endMass, spectraMassArray);
        if(startIdx>0) { startIdx--; }
        if(endIdx<spectraMassArray.length-1) { endIdx++; }

        for(let idx=startIdx; idx<=endIdx; idx++) {
          let s = sortedSpectra[idx];
          let lowerNeutral = (s.pepmass-errorRange) * s.charge - s.charge;
          let upperNeutral = (s.pepmass+errorRange) * s.charge - s.charge;
          if(p.mass<lowerNeutral || p.mass>upperNeutral) { continue; }
          let bFrags = p.fragments['b']['1'];
          let yFrags = p.fragments['y']['1'];

          let bMatches = matchFragments(s.mz, bFrags, config);
          let yMatches = matchFragments(s.mz, yFrags, config);
          
          let bTotal = 0;
          let yTotal = 0;

          for(var k=0; k<bMatches.length; k++) {
            let m = bMatches[k];
            bTotal += s.intensity[m];
          }

          for(var k=0; k<yMatches.length; k++) {
            let m = yMatches[k];
            yTotal += s.intensity[m];
          }

          let totalIntensity = s.intensity.reduce((accumulator, val) => {
            return accumulator + val;
          });

          let score = (bTotal + yTotal) * factorialize(bMatches.length) * factorialize(yMatches.length);
          //let score = bMatches.length + yMatches.length + (bTotal + yTotal)/totalIntensity;

          let spectrumMatch = {
            spectrum: s.title,
            accession: p.accession,
            sequence: p.sequence,
            modifications: p.modifications,
            score: score
          }
          s.matches.push(spectrumMatch);
        }
        
      }
  
    }
    
    process.stdout.write('\n');

    const csvWriter = createCsvWriter({
      path: './old_results.csv',
      header: [
          {id: 'spectrum', title: 'Spectrum title'},
          {id: 'accession', title: 'Accession'},
          {id: 'sequence', title: 'Sequence'},
          {id: 'modifications', title: 'Modifications'},
          {id: 'score', title: 'Score'}
      ]
    });
    let bestMatches = [];
    for(var i=0; i<sortedSpectra.length; i++) {
      let matches = sortedSpectra[i].matches;
      let bestMatch = {};
      bestScore = 0;
      for(var j=0; j<matches.length; j++) {
        if(matches[j].score > bestScore){
          bestScore = matches[j].score;
          bestMatch = matches[j];
        }
      }
      bestMatches.push(bestMatch);
    }

    csvWriter.writeRecords(bestMatches)       // returns a promise
    .then(() => {
        console.log('...Done');
    });
  });
});

function matchFragments(mz, fragments, config) {
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
function closestIdx (num, arr) {
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

function processPeptides(accession, sequence, config) {
  peptidesToReturn = [];
  let cleavedPeptides = cutter.cleave(sequence);
  for(var i=0; i<cleavedPeptides.length; i++) {
    var pep = cleavedPeptides[i];
    let modifications = pepMod.modify(pep.sequence, config.modifications, 2);
    for(var j=0; j<modifications.length; j++) {
      let mass = calculateMass(pep.sequence, modifications[j]);
      if(mass === -1) { continue; }
      if(mass < config.minMass || mass > config.maxMass) { continue; }
      let dbEntry = {
        accession: accession,
        sequence: pep.sequence,
        modifications: modifications[j],
        mass: mass,
        fragments: pepFrag.fragment(pep.sequence, ['b','y'], [1], modifications[j])
      };
      peptidesToReturn.push(dbEntry);
    }
  }
  return peptidesToReturn;
}

function calculateMass(sequence, modifications) {
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
function sortByKey(array, key) {
  return array.sort(function(a, b) {
      var x = a[key]; var y = b[key];
      return ((x < y) ? -1 : ((x > y) ? 1 : 0));
  });
}

// https://medium.freecodecamp.org/how-to-factorialize-a-number-in-javascript-9263c89a4b38
function factorialize(num) {
  if (num < 0) 
        return -1;
  else if (num == 0) 
      return 1;
  else {
      return (num * factorialize(num - 1));
  }
}

// This is from StackOverflow but I'm reading it off my phone so I don't have the full URL. The ID or something is 1374126
Array.prototype.extend = function(other_array) {
  other_array.forEach(function(v) {this.push(v)}, this);
}

function getSum(total, num) {
  return total + num;
}