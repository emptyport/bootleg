'use strict';

let mgf = require('js-mgf');
let fastaParser = require('fasta-js');
let peptideCutter = require('peptide-cutter');
let pepMod = require('peptide-modifier');
let pepFrag = require('peptide-fragmenter');


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
    spectra[i].matches = [];
    let intensitySum = spectra[i].intensity.reduce((total, num) => {
      return total + num;
    });
    for(var j=0; j<spectra[i].intensity.length; j++) {
      spectra[i].intensity[j] = spectra[i].intensity[j] / intensitySum;
    }
  }
  return spectra;
}

const parseFasta = (data) => {
  
}

module.exports.search = (spectra_data, fasta_data, config) => {

  // Let's start off by parsing the spectra and then sorting everything
  let spectra = parseSpectra(spectra_data);
  let sortedSpectra = sortByKey(spectra, 'neutral_mass');
  let spectraLength = sortedSpectra.length;
  let spectraMassArray = new Array(spectraLength);
  for(let i=0; i<spectraLength; i++) {
    spectraMassArray[i] = (sortedSpectra[i].neutral_mass);
  }
  config.minMass = spectraMassArray[0].neutral_mass;
  config.maxMass = spectraMassArray[spectraLength - 1].neutral_mass;

  // Now let's parse the fasta file
  let sequenceDB = fasta.parse(fasta_data);
  let fasta = new fastaParser(config.fastaOptions);
  let cuttingOptions = {
    'enzyme': config.enzyme,
    'num_missed_cleavages': config.missedCleavages,
    'min_length': config.minLength,
    'max_length': config.maxLength
  }
  let cutter = new peptideCutter(cuttingOptions);
  



}

