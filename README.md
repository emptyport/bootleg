# Bootleg

This is a short and sweet implementation of the X! tandem algorithm. A lot of work is still required and things aren't ready for distribution/use yet. Feel free to run `npm install` and take a look at `bootleg.js`. The files are currently hardcoded so you will need to change those. Just provide a fasta file and an mgf file.

As is, Bootleg matches X! tandem about 65% of the time on the test mgf file provided by X! tandem. Not bad, but still a long ways to go. Even though it is single threaded, the runtime on that small test file isn't bad, so I'm optimistic about Bootleg being a viable peptide search engine.