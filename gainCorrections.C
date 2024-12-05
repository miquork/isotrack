#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// Structure to hold the correction factors for each run range
struct RunCorrection {
    int start_run;
    int end_run;
    double corrections[59][7]; // 59 eta indices (including eta=0), 7 depth indices
};

// Function to load the correction data from file
void loadCorrections(const std::string &filename, std::vector<RunCorrection> &run_corrections) {
    std::ifstream file(filename);
    std::string line;

    std::cout << "Reading in gain corrections from file '" << filename << "'\n" << std::flush;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return;
    }

    // Read the header line to parse run ranges
    std::vector<std::pair<int, int>> new_run_ranges;
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token; // Skip the first two headers (#eta dep)
        iss >> token;
        while (iss >> token) {
            int start, end;
            if (std::sscanf(token.c_str(), "%d-%d", &start, &end) == 2) {
                new_run_ranges.emplace_back(start, end);
            }
        }
    }

    // Create RunCorrection objects for each run range
    size_t num_new_run_ranges = new_run_ranges.size();
    std::vector<RunCorrection> new_run_corrections(num_new_run_ranges);
    for (size_t i = 0; i < num_new_run_ranges; ++i) {
        new_run_corrections[i].start_run = new_run_ranges[i].first;
        new_run_corrections[i].end_run = new_run_ranges[i].second;
        // Initialize corrections array with sentinel value
        for (int eta_idx = 0; eta_idx < 59; ++eta_idx) {
            for (int depth_idx = 0; depth_idx < 7; ++depth_idx) {
	      new_run_corrections[i].corrections[eta_idx][depth_idx] = 1.0;//-1.0;
            }
        }
    }

    // Read the data lines and populate corrections
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int eta, depth;
        iss >> eta >> depth;
        std::vector<double> factors;
        double factor;
        while (iss >> factor) {
            factors.push_back(factor);
        }

        // Ensure factors size matches the number of run ranges
        if (factors.size() != num_new_run_ranges) {
            std::cerr << "Mismatch in number of correction factors for eta " << eta << ", depth " << depth << "\n";
            continue;
        }

        // Map eta and depth to indices
        int eta_index = eta + 29;
        int depth_index = depth - 1;

        if (eta_index < 0 || eta_index >= 59 || depth_index < 0 || depth_index >= 7) {
            std::cerr << "Invalid eta or depth value: eta=" << eta << ", depth=" << depth << "\n";
            continue;
        }

        // Assign factors to the corresponding RunCorrection objects
        for (size_t i = 0; i < num_new_run_ranges; ++i) {
            new_run_corrections[i].corrections[eta_index][depth_index] = factors[i];
        }
    }
    file.close();

    // Append new run corrections to the main vector
    run_corrections.insert(run_corrections.end(), new_run_corrections.begin(), new_run_corrections.end());

    // After loading all data, sort the run corrections by start_run
    std::sort(run_corrections.begin(), run_corrections.end(),
        [](const RunCorrection &a, const RunCorrection &b) {
            return a.start_run < b.start_run;
        });
}

// Correction retriever class with caching
class CorrectionRetriever {
public:
 CorrectionRetriever(const std::vector<RunCorrection> &rcs) : run_corrections(rcs), last_index(-1), last_run(-1) {}

    double getCorrection(int run, int eta, int depth) {
        // Corrections only available for [-29,-16] and [16,29] so return 1 elsewhere
        if (abs(eta)>29 || abs(eta)<16) return 1.0;
      
        // Check if run is within the last used run correction
        if (last_index != -1) {
            const RunCorrection &last_rc = run_corrections[last_index];
            if (run >= last_rc.start_run && run <= last_rc.end_run) {
                return retrieveCorrection(last_rc, eta, depth, run);
            }
        }

        // Perform binary search to find the run correction
        auto it = std::lower_bound(run_corrections.begin(), run_corrections.end(), run,
            [](const RunCorrection &rc, int run) {
                return rc.end_run < run;
            });

        if (it == run_corrections.end() || run < it->start_run || run > it->end_run) {
	  if (run!=last_run) {
            std::cerr << "Run " << run << " not found in ranges.\n";
	    if (it != run_corrections.end()) std::cerr << "Nearest range ["<<it->start_run<<","<<it->end_run<<"].\n";
	  }
	  last_run = run;
	  return 1.0;//-1.0;
        }

        last_index = std::distance(run_corrections.begin(), it);
        return retrieveCorrection(*it, eta, depth, run);
    }

private:
    std::vector<RunCorrection> run_corrections; // Store a copy instead of a reference
    int last_index;
    int last_run;

    double retrieveCorrection(const RunCorrection &rc, int eta, int depth, int run) {
        // Map eta and depth to indices
        int eta_index = eta + 29;
        int depth_index = depth - 1;

        if (eta_index < 0 || eta_index >= 59 || depth_index < 0 || depth_index >= 7) {
            std::cerr << "Invalid eta or depth value: eta=" << eta << ", depth=" << depth << "\n";
            return 1.0;//-1.0;
        }

        double correction = rc.corrections[eta_index][depth_index];
        if (correction < 0.0) {
            std::cerr << "No correction factor available for eta " << eta << ", depth " << depth << " in run " << run << ".\n";
            return 1.0;//-1.0;
        }
        return correction;
    }
};

// Global variables
CorrectionRetriever* _gainCorrectionRetriever = nullptr;

int gainCorrections() {
    std::vector<RunCorrection> run_corrections;
    // Load corrections into global run_corrections
    loadCorrections("textfiles/Compensating_Gain_Factors_2024GHI_v5.txt", run_corrections);
    loadCorrections("textfiles/Compensating_Gain_Factors_2024CDEF_New.txt", run_corrections);

    // Initialize global retriever with the global run_corrections
    _gainCorrectionRetriever = new CorrectionRetriever(run_corrections);

    // Example usage
    int run = 383811;
    int eta = -29;
    int depth = 1;

    // Missing runs:
    // 383247

    double correction = _gainCorrectionRetriever->getCorrection(run, eta, depth);

    if (correction > 0.0) {
        std::cout << "Correction factor for run " << run << ", eta " << eta << ", depth " << depth
                  << " is " << correction << "\n";
    }
    std::cout << "Access corrections with _gainCorrectionRetriever->getCorrection(run, ieta, depth)\n";
    
    return 0;
}
