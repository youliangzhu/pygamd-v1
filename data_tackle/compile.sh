g++ -c -std=c++11 MolInfo.cpp XMLParser.cpp XMLBuilder.cpp MSTReader.cpp DCDBuilder.cpp Functions.cpp dataTackle.cpp -I /usr/local/cuda/include/
# nvcc -c -DNVCC rdf.cu strfac.cu 
# nvcc rdf.o strfac.o Functions.o MolInfo.o mst_reader.o Dcdbuilder.o dataTackle.o -o dataTackle
nvcc -c -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -DNVCC rdf.cu strfac.cu
nvcc -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 rdf.o strfac.o Functions.o MolInfo.o XMLParser.o XMLBuilder.o MSTReader.o DCDBuilder.o dataTackle.o -o dataTackle
rm *.o
