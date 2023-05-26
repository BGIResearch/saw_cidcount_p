/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sys/stat.h>
#include <string>
#include <unistd.h>
#include <assert.h>
#include <omp.h>
#include <mutex>
#include <math.h>
#include <iomanip>
#include <map>
#include <vector>
#include <string.h>
#include "hdf5.h"
#include "libdeflate.h"
using namespace::std;
#define DATASETNAME "bpMatrix_"
#define MAJOR_VERSION 1
#define MINOR_VERSION 1
mutex logLock;
void printUsg(string progName);
void printVersionAndAuthor();
// calculate barcode number from file whose format contains .h5 and .bin
uint64_t readOldH5(string maskFile){
  std::string datasetName = DATASETNAME + std::to_string(1);
  hid_t fileID = H5Fopen(maskFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t datasetID = H5Dopen2(fileID, datasetName.c_str(), H5P_DEFAULT);
  
  //read attribute of the dataset
  
  //uint32_t attributeValues[ATTRIBUTEDIM];
  //hid_t attributeID = H5Aopen_by_name(fileID, datasetName.c_str(), ATTRIBUTENAME, H5P_DEFAULT, H5P_DEFAULT);
  //status = H5Aread(attributeID, H5T_NATIVE_UINT32, &attributeValues[0]);
  //cout << "attribute values: " << attributeValues[0] << " "<< attributeValues[1] << endl;
  hid_t dspaceID = H5Dget_space(datasetID);
//  hid_t dtype_id = H5Dget_type(datasetID);
//  hid_t plistID = H5Dget_create_plist(datasetID);
  int rank = H5Sget_simple_extent_ndims(dspaceID);
  hsize_t dims[rank];
  herr_t status = H5Sget_simple_extent_dims(dspaceID, dims, NULL);
  uint64_t bcNum=0;
  uint64_t matrixLen = 1;
  for (int i = 0 ; i<rank; i++){
	matrixLen *= dims[i];
  }
  
  int segment = 1;
  if (rank>=3){
	segment = dims[2];
  }
//  pid_t pd=getpid();
//  cout<<"current memory use "<<get_proc_virtualmem(pd)<<"\t"<<get_proc_virtualmemPeak(pd)<<endl;
  //cout<<"use hyperslab"<<endl;
  uint64_t cpRow=dims[0]/100==0?1:dims[0]/100;
  uint64_t* bpMatrix_buffer2 = new uint64_t[cpRow*dims[1]*dims[2]];
	for(uint64_t rowIdx=0;rowIdx<=dims[0];rowIdx+=cpRow){
	  uint64_t realCpRows=rowIdx+cpRow>dims[0]?dims[0]-rowIdx:cpRow;
	  uint64_t cpSize=realCpRows*dims[1]*dims[2];
	  
	  hsize_t dimsm2[1];
	  dimsm2[0]=cpSize;
	  hid_t memspace2=H5Screate_simple(1,dimsm2,NULL);
	  hsize_t offset[3],count[3];
	  offset[0] = rowIdx;
	  offset[1] = 0;
	  offset[2] = 0;
	  count[0]  = realCpRows;
	  count[1]  = dims[1];
	  count[2]  = dims[2];
	  hsize_t offset_out2[1];
	  hsize_t count_out2[1];
	  offset_out2[0]=0;
	  count_out2[0]=cpSize;
	  if(count[0]*count[1]*count[2]!=cpSize){
		cerr<<"different size\t"<<rowIdx<<"\t"<<count[0]<<"\t"<<count[1]<<"\t"<<cpSize<<endl;
		exit(1);
	  }
	  status = H5Sselect_hyperslab(memspace2, H5S_SELECT_SET, offset_out2, NULL,
								   count_out2, NULL);
	  status = H5Sselect_hyperslab(dspaceID, H5S_SELECT_SET, offset, NULL,
								   count, NULL);
	  status = H5Dread(datasetID, H5T_NATIVE_ULONG, memspace2, dspaceID,
					   H5P_DEFAULT, bpMatrix_buffer2);
	  for (uint32_t r = 0; r < realCpRows; r++){
		//bpMatrix[r] = new uint64_t*[dims[1]];
		for (uint32_t c = 0; c< dims[1]; c++){
		  //bpMatrix[r][c] = bpMatrix_buffer + r*dims[1]*dims[2] + c*dims[2];
		  if (rank >= 3 ){
			segment = dims[2];
			for (int s = 0; s<segment; s++){
			  uint64_t barcodeInt = bpMatrix_buffer2[r*dims[1]*segment + c*segment + s];
			  if (barcodeInt == 0){
				continue;
			  }
			  bcNum++;
			}
		  }else{
			uint64_t barcodeInt = bpMatrix_buffer2[r*dims[1]+c];
			if (barcodeInt == 0){
			  continue;
			}
			bcNum++;
		  }
		}
	  }
	  cout<<rowIdx<<"\t"<<dims[0]<<endl;
	}
  delete[] bpMatrix_buffer2;
	return bcNum;
}
uint64_t readH5(string maskFile){
  herr_t status;
  //open dataset with datasetName
  std::string datasetName = DATASETNAME + std::to_string(1);
  hid_t fileID = H5Fopen(maskFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t datasetID = H5Dopen2(fileID, datasetName.c_str(), H5P_DEFAULT);
  if(datasetID==H5I_INVALID_HID){
    cerr<<"Error, fail to open h5 file"<<endl;
    exit(1);
  }
//#ifdef OPT_MEMORY_PLAN_B
//  cout<<"before reserve:\t"<<bpMap.getAllocatedMemorySize()<<endl;
//  uint64_t bcNum=5073598806;
//  bpMap.reserve(bcNum);
//  cout<<"after reserve:\t"<<bpMap.getAllocatedMemorySize()<<endl;
//#endif
  hid_t dspaceID = H5Dget_space(datasetID);

//  hid_t dtype_id = H5Dget_type(datasetID);
  
  hid_t plistID = H5Dget_create_plist(datasetID);
  
  int rank = H5Sget_simple_extent_ndims(dspaceID);
  
  hsize_t dims[rank];
  status = H5Sget_simple_extent_dims(dspaceID, dims, NULL);
  
  
  uint64_t matrixLen = 1;
  for (int i = 0; i < rank; i++) {
	matrixLen *= dims[i];
  }
  
  int segment = 1;
  if (rank >= 3) {
	segment = dims[2];
  }
  
  status = H5Sselect_all(dspaceID);
  hsize_t nchunks;
  status=H5Dget_num_chunks(datasetID, dspaceID, &nchunks);
  if(status<0){
	return status;
  }
  size_t chunk_dims[rank];
  int flag = 0;
  for (int r = 0; r < rank; r++) {
	if (dims[r] == 1) {
	  flag++;
	  chunk_dims[r] = 1;
	} else {
	  chunk_dims[r] = 0;
	}
	
  }
  hsize_t offset_next[rank];
  for (uint64_t cid = 1; cid < nchunks; cid++) {
	H5Dget_chunk_info(datasetID, dspaceID, cid, offset_next, NULL, NULL, NULL);
//	cout<<cid<<"\t"<<offset_next[0]<<"\t"<<offset_next[1]<<"\t"<<offset_next[2]<<endl;
	for (int r = 0; r < rank; r++) {
	  if (!chunk_dims[r] && offset_next[r]) {
		chunk_dims[r] = offset_next[r];
		flag++;
	  }
	}
	if (flag >= rank) {
	  break;
	}
  }
  for (int r = 0; r < rank; r++){
    if(chunk_dims[r] == 0){
	  chunk_dims[r] = dims[r];
    }
  }
  hsize_t chunk_len = 1;
  for (int r = 0; r < rank; r++) {
	chunk_len *= chunk_dims[r];
  }

  hsize_t chunkBufSize=10;
  hsize_t **compressed_buffer = new hsize_t *[chunkBufSize];
  hsize_t **buffer = new hsize_t *[chunkBufSize];
  hsize_t **offset = new hsize_t *[chunkBufSize];
  hsize_t **compressed_buffer2 = new hsize_t *[chunkBufSize];
  hsize_t **buffer2 = new hsize_t *[chunkBufSize];
  hsize_t **offset2 = new hsize_t *[chunkBufSize];
  for (hsize_t i = 0; i < chunkBufSize; i++) {
	compressed_buffer[i] = new hsize_t[chunk_len];
	buffer[i] = new hsize_t[chunk_len];
	offset[i] = new hsize_t[rank];
	compressed_buffer2[i] = new hsize_t[chunk_len];
	buffer2[i] = new hsize_t[chunk_len];
	offset2[i] = new hsize_t[rank];
  }
  hsize_t *chunk_size = new hsize_t[nchunks];
  
  uint64_t bcNum=0;
  uint64_t chunk_num = 0;
  uint64_t now_chunk[3];
  bool is_complete_hdf5 = false;
  string buf1status="empty";
  string buf2status="empty";
  uint64_t buf1size=0;
  uint64_t buf2size=0;
  int openMPthreads=3;
  openMPthreads=2;
//  cout<<"total chunks number:\t"<<nchunks<<endl;
#pragma omp parallel num_threads(openMPthreads)
  {
	int num_id = omp_get_thread_num();
	if (num_id == 0) {
	  libdeflate_decompressor *decompressor = libdeflate_alloc_decompressor();
	  now_chunk[0]=0;
	  uint64_t chunk_start=0;
	  uint64_t chunk_end=chunk_start+chunkBufSize>nchunks?nchunks:chunk_start+chunkBufSize;
//	  cout<<"0\t"<<now_chunk[0]<<"\t"<<nchunks<<endl;
	  while(now_chunk[0]<nchunks){
	    if(buf1status=="empty"){
//	      cout<<chunk_start<<"\t"<<chunk_end<<endl;
		  for(uint64_t chunk_index=chunk_start;chunk_index<chunk_end;chunk_index++){
		    uint64_t patchIdx=chunk_index-chunk_start;
			uint32_t filter = 0;
			size_t actual_out = 0;
			if(H5Dget_chunk_info(datasetID, dspaceID, chunk_index, offset[patchIdx], &filter, NULL,&chunk_size[patchIdx])<0){
			  cerr<<"Error, get chunk info error, please check h5 file"<<endl;
			  exit(1);
			}
			if(H5Dread_chunk(datasetID, H5P_DEFAULT, offset[patchIdx], &filter,compressed_buffer[patchIdx])<0){
			  cerr<<"Error, read chunk info error, please check h5 file"<<endl;
			  exit(1);
			}
			int dstatus = libdeflate_zlib_decompress(decompressor, (void *) compressed_buffer[patchIdx],chunk_size[patchIdx], (void *) buffer[patchIdx],chunk_len * sizeof(uint64_t), &actual_out);
			chunk_num++;
			now_chunk[0]++;
		  }
		  buf1size=chunk_end-chunk_start;
		  assert(buf1size<=chunkBufSize);
		  buf1status="full";
		  chunk_start=chunk_end;
		  chunk_end=chunk_start+chunkBufSize>nchunks?nchunks:chunk_start+chunkBufSize;
		}
		if(buf2status=="empty"){
		  if(buf1status=="full"){
		    hsize_t **tmpPtr1=compressed_buffer2;
		    hsize_t **tmpPtr2=buffer2;
		    hsize_t **tmpPtr3=offset2;
		    compressed_buffer2=compressed_buffer;
		    buffer2=buffer;
		    offset2=offset;
		    compressed_buffer=tmpPtr1;
		    buffer=tmpPtr2;
		    offset=tmpPtr3;
		    buf1status="empty";
		    buf2status="full";
		    buf2size=buf1size;
		    assert(buf2size<=chunkBufSize);
		    buf1size=0;
		  }
		}else{
		  usleep(10);
		}
		
	  }
	  while(1){
		if(buf2status=="empty"){
		  if(buf1status=="full"){
		    hsize_t **tmpPtr1=compressed_buffer2;
		    hsize_t **tmpPtr2=buffer2;
		    hsize_t **tmpPtr3=offset2;
		    compressed_buffer2=compressed_buffer;
		    buffer2=buffer;
		    offset2=offset;
		    compressed_buffer=tmpPtr1;
		    buffer=tmpPtr2;
		    offset=tmpPtr3;
		    buf1status="empty";
		    buf2status="full";
		    buf2size=buf1size;
		    assert(buf2size<=chunkBufSize);
		    buf1size=0;
		  }else{
		    break;
		  }
		  usleep(10);
		}
	  }
//	  logLock.lock();
//	  cout<<"0 thread complete\t"<<now_chunk[0]<<"\t"<<nchunks<<endl;
//	  logLock.unlock();
	  is_complete_hdf5 = true;
	  libdeflate_free_decompressor(decompressor);
	  
	  status = H5Dclose(datasetID);
	  status = H5Fclose(fileID);
	}
	if (num_id == 1) {
	  // HashTable
	  now_chunk[1] = 0;
	  uint64_t lastValue=0;
	  while(now_chunk[1] < nchunks){
		if(buf2status=="full") {
		  for (uint64_t i = 0; i < buf2size; i++) {
			for (uint32_t y = offset2[i][0]; y < min(offset2[i][0] + chunk_dims[0], dims[0]); y++) {
			  for (uint32_t x = offset2[i][1]; x < min(offset2[i][1] + chunk_dims[1], dims[1]); x++) {
			    if(rank>=3) {
				  for (uint32_t z = offset2[i][2]; z < min(offset2[i][2] + chunk_dims[2], dims[2]); z++) {
					uint64_t barcodeInt = buffer2[i][(y - offset2[i][0])*chunk_dims[1]*dims[2] + (x - offset2[i][1])*dims[2] + z];
					if(barcodeInt==0) {
					  continue;
					}
					bcNum++;
				  }
				}else{
				  uint64_t barcodeInt = buffer2[i][(y - offset2[i][0])*chunk_dims[1] + x - offset2[i][1]];
				  if(barcodeInt==0) {
					continue;
				  }
				  bcNum++;
			    }
			  }
			}
			now_chunk[1]++;
		  }
//		  logLock.lock();
//		  cout<<"1 thread processed chunks:\t"<<now_chunk[1]<<". bcNum:\t"<<bcNum<<endl;
//		  logLock.unlock();
		  buf2status = "empty";
		} else {
		  usleep(10);
		}
	  }
//	  logLock.lock();
//	  cout<<"1 thread complete\t"<<now_chunk[1]<<"\t"<<nchunks<<endl;
//	  logLock.unlock();
	}
  }
  uint64_t chunkSize=0;
  chunkSize=chunkBufSize;
  for (uint32_t i = 0; i < chunkSize; i++) {
	delete[] compressed_buffer[i];
	delete[] buffer[i];
	delete[] offset[i];
	delete[] compressed_buffer2[i];
	delete[] buffer2[i];
	delete[] offset2[i];
  }
  delete[] chunk_size;
  delete[] compressed_buffer;
  delete[] buffer;
  delete[] offset;
  delete[] compressed_buffer2;
  delete[] buffer2;
  delete[] offset2;
  return bcNum;
}
int main(int argc,char* argv[]) {
	// -i <h5 file>
	// -t <threads number|optional>
	// -o <output file|optional>
  if(argc<4) {
    if(argc==2) {
	  if(strcmp(argv[1], "-h")!=0) {
		cerr << "Arguments error!" << endl;
	  }
	}else{
	  cerr << "Arguments error!" << endl;
    }
	printUsg(string(argv[0]));
  }
  // process parameters
  string maskFile="";
  int threadsNum=16;
  string outputFile="";
  uint64_t bcNum=0;
  string species="";
  string genomeMemoryFile="";
  double genomeFileSize=0;
  const char *shortOpt="i:s:r:g:o:h";
  int nextOpt;
  while(-1!=(nextOpt=getopt(argc,argv,shortOpt))){
    switch(nextOpt){
      case 'i':maskFile=optarg;break;
      case 's':species=optarg;break;
      case 'r':genomeMemoryFile=optarg;break;
      case 'g':genomeFileSize=atof(optarg);break;
//      case 't':threadsNum=atoi(optarg);break;
      case 'o':outputFile=optarg;break;
//      case 'h':printUsg(string(argv[0]));break;
      default:{
        cerr<<"Error, no such parameter,"<<char(nextOpt)<<endl;
        exit(1);
      }
    }
  }
  // check input file
  ifstream maskIn(maskFile);
  if(!maskIn){
    cerr<<"Error, cannot open such file,"<<maskFile<<endl;
    exit(1);
  }
  maskIn.close();
  
  
  
  if(maskFile.rfind(".bin")==maskFile.size()-4){
    // bin
    struct stat fileInfo;
    if(stat(maskFile.c_str(),&fileInfo)==0){
    	uint64_t fileSize=fileInfo.st_size;
    	if(fileSize%16!=0){
    	  cerr<<"Error, mask file was not complete, check the file size("<<fileSize<<" bytes)"<<endl;
    	  exit(1);
    	}
	  	bcNum=fileSize/16;
    }else{
      cerr<<"Error, fail to stat the file"<<endl;
      exit(1);
    }
  }else if(maskFile.rfind(".h5")==maskFile.size()-3){
    // h5
    // single thread
	bcNum=readH5(maskFile);
//	bcNum=readOldH5(maskFile);
  }else{
    cerr<<"Error, mask file should be in .bin or .h5 format"<<endl;
    exit(1);
  }
//  for(uint64_t i=100*1024*1024;i<=1*1024*1024*1024;i+=100*1024*1024){
//	uint64_t needMemory=1ul<<uint64_t(ceil(log2((double)i/12))+8);
//    cout<<i<<"\t"<<ceil(log2(i/12))+8<<"\t"<<needMemory<<"\t"<<needMemory/(1024*1024*1024)<<endl;
//  }
//  for(uint64_t i=402653184;i<402653186;i++){
//	uint64_t needMemory=1ul<<uint64_t(ceil(log2((double)i/12))+8);
//	cout<<i<<"\t"<<ceil(log2(i/12))+8<<"\t"<<needMemory<<"\t"<<needMemory/(1024*1024*1024)<<endl;
//  }
//  for(uint64_t i=805306367;i<=805306369;i++){
//	uint64_t needMemory=1ul<<uint64_t(ceil(log2((double)i/12))+8);
//	cout<<i<<"\t"<<ceil(log2(i/12))+8<<"\t"<<needMemory<<"\t"<<needMemory/(1024*1024*1024)<<endl;
//  }
  uint64_t bcMemory=1ul<<(uint64_t(ceil(log2(bcNum/12)))+8);
  bcMemory=(double)bcMemory/(1024*1024*1024);
  if(bcMemory<2){
	bcMemory=1;
  }
  double refMemory=0;
  map<string,double > ref2mem;
  if(genomeMemoryFile!="") {
	ifstream refIn(genomeMemoryFile);
	if(!refIn) {
	  cerr << "Error, cannot open such file," << genomeMemoryFile << endl;
	  exit(1);
	}
	string tmpLine;
	while(getline(refIn,tmpLine)){
	  vector<string> eles;
	  string tmpEle="";
	  for(uint64_t i=0;i<tmpLine.size();i++){
		if(tmpLine[i]=='\t'){
		  eles.emplace_back(tmpEle);
		  tmpEle="";
		}else{
		  tmpEle+=tmpLine[i];
		}
	  }
	  if(!tmpEle.empty()){
		eles.emplace_back(tmpEle);
	  }
	  if(eles.size()!=2){
		cerr<<"Error, format of genome memory file error which should contain 2 columns"<<endl;
		exit(1);
	  }
	  ref2mem[eles[0]]=atof(eles[1].c_str());
	}
	refIn.close();
  }
  
  if(ref2mem.find(species)!=ref2mem.end()){
    refMemory=ceil(ref2mem[species]);
  }else{
    refMemory=genomeFileSize*12;
  }
  uint64_t totalMemeory=bcMemory+ceil(refMemory)+10;
  if(outputFile.empty()){
    cout<<bcNum<<endl;
    cout<<totalMemeory<<endl;
  }else{
	ofstream out(outputFile);
    out<<bcNum<<endl;
    out<<totalMemeory<<endl;
	out.close();
  }
  return 0;
}
void printUsg(string progName) {
//  cerr<<"Arguments error!"<<endl;
printVersionAndAuthor();
  cerr<<progName<<"  [options]"<<endl;
  cerr<<"\t"<<"-i <h5 file>"<<endl;
  cerr<<"\t"<<"-s <species>"<<endl;
  cerr<<"\t"<<"-g <Genome file size of the species in STAR index(GB)>"<<endl;
  cerr<<"\t"<<"-r <genome memory file|optional>"<<endl;
//	cerr<<"\t"<<"-t <threads number|optional>"<<endl;
  cerr<<"\t"<<"-o <output file|optional>, default output to stdout"<<endl;
  cerr<<"\t"<<"-h print help"<<endl;
  exit(1);
}
void printVersionAndAuthor(){
  cout<<"Program: calBcNum\n";
  cout<<"Version: "<<MAJOR_VERSION<<"."<<MINOR_VERSION<<endl;
  cout<<"Contact: GongChun<gongchun@genomics.cn>"<<endl;
}
