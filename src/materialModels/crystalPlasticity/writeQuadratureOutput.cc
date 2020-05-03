#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::writeQuadratureOutput(std::string _outputDirectory, unsigned int _currentIncrement)
{
	this->pcout << "writing Quadrature data to file\n";
	//
	//set output directory, if provided
	std::string dir(_outputDirectory);

	if (_outputDirectory.back() != '/') dir += "/";

	std::string fileName("QuadratureOutputs");
	std::string fileExtension(".csv");
	fileName += std::to_string(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
	std::ofstream file((dir + fileName + std::to_string(_currentIncrement)+fileExtension).c_str());
	char buffer[200];
	if (file.is_open()) {
		for (std::vector<std::vector<double> >::iterator it = outputQuadrature.begin(); it != outputQuadrature.end(); ++it) {
			for (std::vector<double>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				sprintf(buffer, "%8.5e ,", *it2);
				file << buffer;
			}
			file << std::endl;
		}
		file.close();
	}
	else {
		this->pcout << "Unable to open file for writing quadrature outputs\n";
		exit(1);
	}

	//join files from all processors into a single file on processor 0
	//and delete individual processor files
	MPI_Barrier(MPI_COMM_WORLD);
	if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
		std::string fileName2("QuadratureOutputs");
		std::ofstream file2((dir + fileName2 + std::to_string(_currentIncrement)+fileExtension).c_str());
		for (unsigned int proc = 0; proc<dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); proc++) {
			std::string fileName3("QuadratureOutputs");
			fileName3 += std::to_string(proc);
			std::ofstream file3((dir  + fileName3 + std::to_string(_currentIncrement)+fileExtension).c_str(), std::ofstream::in);
			file2 << file3.rdbuf();
			//delete file from processor proc
			remove((dir + fileName3 + std::to_string(_currentIncrement)+fileExtension).c_str());
		}
		file2.close();
	}
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
