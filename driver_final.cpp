#include "PolyRing.h"

using namespace std;

int main(int argc, char *argv[]) {
	
	if (argc > 3) {
		cout << "Too many parameter. Exit";
		return 0;
	}
	
	int specified_param = -1;
	if (argc == 2) {
		specified_param = atoi(argv[1]);
	}
	
	vector<pair<int, int>> security_params = {{128, 200000}, {192, 350000}, {256, 512000}};
	
	int message_space_size = 256;
	int run_amount = 10000;
	int seed = 0;
	
	for (auto check : security_params) {
		int S = check.first;
		int high_n = check.second;
		
		if (specified_param != -1 && specified_param != S) {
			continue;
		}
		
		string freq_file_name = "freq_output" + to_string(S) + ".csv";
		string res_file_name = "res_output" + to_string(S) + ".csv";
		ofstream res_file(res_file_name.c_str(), ofstream::out); 
		cout << "Checking for security param " << S << endl;
		
		for(int l = 1; l <= 256; l *= 2) {
			int h = (S + l - 1) / l;
			int k = (message_space_size + l - 1) / l;
			
			
			int low_n = 0;
			while (low_n <= high_n - 16) {
				int n = (low_n + high_n) / 2;
				cout << "Running for l = " << l << ", n = " << n << ", h = " << h << ", k = " << k << endl; 
				
				MersennePoly C(n, k, h, l, seed);
				auto result = C.runs_final(run_amount, freq_file_name.c_str());
				
				double accuracy = result[0];
				double avg_ms = result[1];
				double erfc1 = result[2];
				double erfc2 = result[3];
				double avg_flips = result[4];
				
				int key_size = 2 * n * l;
				double key_ratio = double(1513678) / double(key_size);
				double security_value = log2((3.0 * l) / (4 * n) + 2.0 / 3) * 3 * h * l;
				
				int pass = int(accuracy * run_amount);

				res_file << n << ", " << k << ", " << h <<", " << l << ", " << key_size << ", ";
				res_file << key_ratio << ", " << erfc1 << ", " << erfc2 << ", ";
				res_file << run_amount << ", " << pass << ", " << accuracy << ", " << avg_ms << endl;
				
				
				if (accuracy >= 0.99) {
					high_n = n - 1;
				} else {
					low_n = n + 1;
				}
			}
			high_n = high_n * 2;
		}
		
		res_file.close();
	}
    return 0;
}
