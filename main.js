#!/usr/bin/env node
'use strict';

var pjson = require('./package.json');
const minimist = require('minimist');
let fs = require('fs');
let bootleg = require('./bootleg');

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
  console.log('  -f\tspecify the fasta file (in .fasta format) - required');
  console.log('  -s\tspecify the spectra file (in .mgf format) - required');
  console.log('  -c\tspecify the config file (in .json format)\n\t- defaults to default_config.json');
  console.log('  -o\tspecify the output file (in .csv format)\n\t- defaults to the name of the spectra file');
  console.log("\nYou are using Bootleg v"+pjson.version);
}

console.log('\n<===Running Bootleg===>');
console.log('Reading in fasta file...');
fs.readFile(args.f, 'utf-8', function(fasta_err, fasta_data) {
  if (fasta_err) throw fasta_err;
  
  console.log('Reading in spectra file...');
  fs.readFile(args.s, 'utf-8', function(spectra_err, spectra_data) {
    if (spectra_err) throw spectra_err;
    
    console.log('Reading in config...');
    let config = require(args.c);

    console.log('Running search...');
    bootleg.search(spectra_data, fasta_data, config);
    
    console.log('Finished!');
  });
});