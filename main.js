#!/usr/bin/env node
'use strict';

var pjson = require('./package.json');
const minimist = require('minimist');
let fs = require('fs');
let bootleg = require('./bootleg');
const createCsvWriter = require('csv-writer').createObjectCsvWriter;

let args = minimist(process.argv.slice(2), {
  default: {
    c: './default_config.json'
  },
  alias: {
    h: 'help',
    v: 'version'
  }
});

if(args.v) {
  console.log("You are using Bootleg v"+pjson.version);
}

if(args.h) {
  console.log("Usage: node main.js [options]\n");
  console.log("Options");
  console.log('  -s\tspecify the spectra file (in .mgf format) - required');
  console.log('  -c\tspecify the config file (in .json format)\n\t- defaults to default_config.json');
  console.log('  -o\tspecify the output file (in .csv format)\n\t- defaults to the name of the spectra file');
  console.log("\nYou are using Bootleg v"+pjson.version);
}

console.log('\n<===Running Bootleg===>');

console.log('Reading in spectra file...');
fs.readFile(args.s, 'utf-8', function(spectra_err, spectra_data) {
  if (spectra_err) throw spectra_err;
  
  console.log('Reading in config...');
  let config = require(args.c);

  console.log('Running...');
  let matches = bootleg.search(spectra_data, config);

  console.log('Saving results...');
  const csvWriter = createCsvWriter({
    path: args.o,
    header: [
        {id: 'spectrum', title: 'Spectrum title'},
        {id: 'accession', title: 'Accession'},
        {id: 'sequence', title: 'Sequence'},
        {id: 'modifications', title: 'Modifications'},
        {id: 'score', title: 'Score'}
    ]
  });

  let formattedResults = [];
  for(let i=0; i<matches.matches.length; i++) {
    let match = matches.matches[i];
    formattedResults.push({
      spectrum: match.spectrum,
      accession: match.accession,
      sequence: match.sequence,
      modifications: match.modifications,
      score: match.bestScore
    });
  }
  csvWriter.writeRecords(formattedResults).then(() => {
    console.log('Finished!');
  });

});
