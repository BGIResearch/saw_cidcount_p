# calBcNum
calculate CID number in mask file and evaluate an exact memory usage in SAW_mapping
## Install
```shell
git clone git@github.com:BGIResearch/SAW_CIDcount.git
cd SAW_CIDcount
make
```
## Run
```text
Program: calBcNum
Version: 1.1
Contact: GongChun<gongchun@genomics.cn>
calBcNum  [options]
	-i <h5 file>
	-s <species>
	-g <Genome file size of the species in STAR index(GB)>
	-r <genome memory file|optional>
	-o <output file|optional>, default output to stdout
	-h print help
```