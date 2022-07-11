Tthis is project for designing a scaffolding tool with miBF usage. The algorithm imitates LINKS as much as it can.

Submodules:
  * btllib
  * btl_bloomfilter

Dependent subdirectory:
  * Common

After cloning directly from the repository run:
```bash
./autogen.sh
```
Compiling BioBloomTools should be as easy as:
```bash
./configure && make
```

To build miBF
```bash
MIBFConstruct/biobloommimaker draft.fa
```

To map miBF
```bash
MIBFAnalyze/mibfanalyzzereads -m /path/mibfdraft.bf -r reads.fa 
```
