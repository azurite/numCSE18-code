#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

VectorXd load_pgm(const std::string &filename) {
	// returns a picture as a flattened vector

	int row = 0, col = 0, rows = 0, cols = 0;

	std::ifstream infile(filename);
	std::stringstream ss;
	std::string inputLine = "";

	// First line : version
	std::getline(infile,inputLine);

	// Second line : comment
	std::getline(infile,inputLine);

	// Continue with a stringstream
	ss << infile.rdbuf();
	// Third line : size
	ss >> cols >> rows;

	VectorXd picture(rows*cols);

	// Following lines : data
	for(row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			int val;
			ss >> val;
			picture(col*rows + row) = val;
		}
	}

	infile.close();
	return picture;
}

// example usage: "a.out 02 sad"
// will identify subject02.sad.pgm
int main(int argc, char **argv) {

	int h = 231;
	int w = 195;
	int M = 15;

	MatrixXd faces(h*w, M);
	VectorXd meanFace = VectorXd::Zero(h*w);

	// loads pictures as flattened vectors into faces
	for (int i=0; i<M; i++) {
		std::string filename = "./basePictures/subject"+
			std::to_string(i+1) + ".pgm";
		VectorXd flatPic = load_pgm(filename);
		faces.col(i) = flatPic;

		meanFace += flatPic;
	}

	meanFace /= M;
	faces.colwise() -= meanFace;

	JacobiSVD<MatrixXd> svd(faces, ComputeThinU | ComputeThinV);
	MatrixXd eigenfaces = svd.matrixU();

	// try to recognize a test face
	string testPicName = "./testPictures/subject" + string(argv[1]) + "." + string(argv[2]) + ".pgm";
	VectorXd newFace = load_pgm(testPicName);

	// TODO: Point (f)
	VectorXd projNewFace = eigenfaces.transpose() * (newFace - meanFace);
	MatrixXd projFaces = eigenfaces.transpose() * faces;

	int k;
	(projFaces.colwise() - projNewFace).colwise().norm().minCoeff(&k);
	std::cout << testPicName << " identified as subject " << k + 1 << std::endl;
}
