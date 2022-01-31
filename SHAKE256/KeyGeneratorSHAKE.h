#include <cryptopp/cryptlib.h>
#include <cryptopp/shake.h>

#include <NTL/tools.h>
#include <NTL/pair.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>

#include <bits/stdc++.h>

using namespace std;
using namespace NTL;
using namespace CryptoPP;

const int seed_len = 32;
struct KeyGeneratorSHAKE{
	SHAKE256 shake;
	int log_N;
	int N, H, L;
	unsigned char seed[seed_len];
	 
	int coeff_byte_len;
	unsigned char* coeff_byte_arr;
	
		
	///restart + update shake with the fixed seed  
	void restart() {
		shake.Restart();
		shake.Update(seed, seed_len);
	}
	
	void setSeed(const unsigned char* seed) {
		copy(seed, seed + seed_len, this->seed); 
	}
	
	KeyGeneratorSHAKE(int N, int H, int L) {
		this->N = N;
		this->H = H;
		this->L = L;
		
		
		
		log_N = 0;
		while ((1LL << log_N) < N) log_N++;
		
		coeff_byte_len = (N + 7) / 8;
		coeff_byte_arr = new unsigned char[coeff_byte_len];
	}
	~KeyGeneratorSHAKE() {
		delete[] coeff_byte_arr;
	}
		
	
	int getRandomBit(int& use_buffer_bit, unsigned char *buffer, int buffer_len) {
		int res = 0;
		int need = log_N;
		while (need > 0) {
			int rem_byte = 8 - (use_buffer_bit % 8);
			int use = min(need, rem_byte);
			res = (res << use) + ((int(buffer[use_buffer_bit / 8]) % (1 << rem_byte)) >> (rem_byte - use));
			use_buffer_bit += use;  
			need -= use; 
		}
		return res;
	}
	
	
	///expect ZZ_p and ZZ_pE initialized correctly before calling these 
	ZZ_p GenerateHSparseP(int &use_buffer_bit, unsigned char *buffer, int buffer_len) {
		memset(coeff_byte_arr, 0, coeff_byte_len);
		for(int i = 1; i <= H; ) {
			int cur_bit = getRandomBit(use_buffer_bit, buffer, buffer_len);
			if (cur_bit >= N) continue;
			
			int byte_pos = cur_bit / 8;
			int bit_pos = cur_bit % 8;
			if (!(coeff_byte_arr[byte_pos] & (1 << bit_pos))) {
				coeff_byte_arr[byte_pos] |= (1 << bit_pos);
				i++;
			}
		}
		ZZ coeff = ZZFromBytes(coeff_byte_arr, coeff_byte_len);
		return conv<ZZ_p>(coeff);
	}
	 
	ZZ_pE GenerateHSparse() {
		int copy = 16;
		int avg_bytes = ((log_N * H * L + 7) / 8);
		
		while (1) {
			//assert(copy <= 1024);			
			int buffer_len = avg_bytes * (copy + 2); /// by Markov, <= 1 / copy chance that this will fail 
			unsigned char* buffer = new unsigned char[buffer_len];
			restart();
			shake.TruncatedFinal(buffer, buffer_len);
			
			
			ZZ_pX P;
			
			int use_buffer_bit = 0;
			bool usable = true;
			for(int i = 0; i < L; i++) {
				if (use_buffer_bit > avg_bytes * 8 * copy) {
					usable = false;
					break;
				}
				SetCoeff(P, i, GenerateHSparseP(use_buffer_bit, buffer, buffer_len));
				
			}
			
			delete[] buffer;
			if (usable) {
				return conv<ZZ_pE>(P);
			}
		}
	}
	
	ZZ_pE GenerateRandom() {
		int required_bit = N * L;
		int required_bytes = (required_bit + 7) / 8;
		//cout << "Alloc " << required_bytes << endl;
		unsigned char* buffer = new unsigned char[required_bytes];
		restart();
		shake.TruncatedFinal(buffer, required_bytes);
		
		
		
		
		ZZ_pX res; 
		int use_bit = 0;
		for(int i = 0; i < L; i++) {
			///realign the bits 
			memset(coeff_byte_arr, 0, coeff_byte_len);
			int cur_bit = 0; 
			
			while (cur_bit < N) {
				int rem_use_byte = 8 - (use_bit % 8);
				int rem_cur_byte = 8 - (cur_bit % 8); 
				int use = min(min(rem_use_byte, rem_cur_byte), N - cur_bit);
				
				
				unsigned char pull_val = ((int(buffer[use_bit / 8]) % (1 << rem_use_byte)) >> (rem_use_byte - use));
				coeff_byte_arr[cur_bit / 8] |= pull_val << (8 - rem_cur_byte);  
				use_bit += use;  
				cur_bit += use; 
			}
			//assert(cur_bit == N);
			ZZ coeff = ZZFromBytes(coeff_byte_arr, coeff_byte_len);
			SetCoeff(res, i, conv<ZZ_p>(coeff));
		}
		
		
		delete[] buffer;
		return conv<ZZ_pE>(res);
	}
};

