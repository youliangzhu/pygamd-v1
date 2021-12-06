g++ -c -std=c++11 MolInfo.cpp mst_reader.cpp Dcdbuilder.cpp Functions.cpp dataTackle.cpp -I /usr/local/cuda/include/ 
nvcc -c -DNVCC rdf.cu strfac.cu 
nvcc rdf.o strfac.o Functions.o MolInfo.o mst_reader.o Dcdbuilder.o dataTackle.o -o galaTackle
rm *.o
