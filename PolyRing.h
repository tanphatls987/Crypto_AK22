#include <NTL/tools.h>
#include <NTL/pair.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>

#include <chrono>
#include <fstream>
#include <vector>

#include "SHAKE256/KeyGeneratorSHAKE.h"

using namespace std;
using namespace NTL;

class MersennePoly {
    int N, K, H, L;
    int rand_seed;
    int flips;
    mt19937_64 mersenne_twister;
    KeyGeneratorSHAKE key_generator;
    
    int RHO;
    double *freq[2];

	unsigned char *HP_bytes; /// used by GenerateHSparseP
	unsigned char *message;
	unsigned char **encoded_message; /// repeating every bit RHO times
	unsigned char *recovered_message;
	ZZ_pE c1;
	unsigned char **c2;
	unsigned char **decoder_inp; 
	
    ZZ_p GenerateHSparseP() { 
        memset(HP_bytes, 0, N/8 + 1);
        int i;
        for(i = 0; i < H/8; i++) {
            HP_bytes[i] = 0xFF;
        }
        HP_bytes[i] = (1U << (H % 8)) - 1;
        for(i = H - 1; i >= 0; i--) {
            int pos = RandomBnd(N - i);
            int a = i, b = i + pos;
            int val1 = (HP_bytes[a/8]>>(a%8))&1;
            int val2 = (HP_bytes[b/8]>>(b%8))&1;
            int xor_val = val1 ^ val2;
            HP_bytes[a/8] ^= xor_val<<(a%8);
            HP_bytes[b/8] ^= xor_val<<(b%8);
        }
        ZZ_p res = conv<ZZ_p>(ZZFromBytes(HP_bytes, N/8 + 1));
        //cout << "Weight of H Sparse --> "<< weight(conv<ZZ>(x)) << endl;
        return res;
    }

	
    ZZ_pE GenerateRandom() {
        ZZ_pX P;
        for(int i = 0; i < L; i++) {
            SetCoeff(P, i, random_ZZ_p());
        }
        return conv<ZZ_pE>(P);
    }

    ZZ_pE GenerateHSparse() {
        ZZ_pX P;
        for(int i = 0; i < L; i++) {
            SetCoeff(P, i, GenerateHSparseP());
        }
        return conv<ZZ_pE>(P);
    }
    
    ZZ_pE GenerateRandom_SHAKE() {
        long long shake_seed[8]; /// 64 * 4 = 256 bits
		for(int i = 0; i < 8; i++) shake_seed[i] = mersenne_twister();
		key_generator.setSeed((unsigned char*)shake_seed);
		return key_generator.GenerateRandom();
    }

    ZZ_pE GenerateHSparse_SHAKE() {
		long long shake_seed[8]; /// 64 * 4 = 256 bits
		for(int i = 0; i < 8; i++) shake_seed[i] = mersenne_twister();
		key_generator.setSeed((unsigned char*)shake_seed);
		return key_generator.GenerateHSparse();
    }
    
    

    void MessageEncoderP(unsigned char* bytes, unsigned char *message, int bit_position) {
        // input -> 1010, RHO = 4
        // output -> 1111000011110000
        for(int i = 0; i < K; i++) {
            for(int j = 0; j < RHO; j++) {
                int pos1 = i;
                int pos2 = i * RHO + j;
                int val1 = (message[(bit_position + pos1)/8]>>((bit_position + pos1)%8))&1;
                bytes[pos2/8] ^= val1<<(pos2%8);
            }
        }
        //ZZ_p res = conv<ZZ_p>(ZZFromBytes(bytes, N/8 + 1));
    }

    void MessageEncoder(unsigned char **bytes, unsigned char *message) {
        //ZZ_pX P;
        for(int i = 0; i < L; i++) {
            memset(bytes[i], 0, N/8 + 1);
            MessageEncoderP(bytes[i], message, K*i);
            //SetCoeff(P, i, MessageEncoderP(message, K*i));
        }
    }

    void XORP(unsigned char* bytes, ZZ_p a, unsigned char *bbB) {
        ZZ aa = conv<ZZ>(a);
        BytesFromZZ(bytes, aa, N/8 + 1);
        for(int i = 0; i < N/8 + 1; i++) {
            bytes[i] ^= bbB[i];
        }
    }

    void XOR(unsigned char **bytes, ZZ_pE a, unsigned char** b) {
        ZZ_pX P;
        ZZ_pX aa = conv<ZZ_pX>(a);
        for(int i = 0; i < L; i++) {
            XORP(bytes[i], coeff(aa, i), b[i]);
        }
    }

    void MessageDecoderP(unsigned char *plainText, int bit_position, unsigned char* bytes, unsigned char *original_message) {
        for(int i = 0; i < K; i++) {
            int ones = 0;
            for(int j = i*RHO; j < (i + 1)*RHO; j++) {
                int val = (bytes[j/8]>>(j%8))&1;
                ones += val;
            }
            if(ones > RHO/2) {
                plainText[(bit_position + i)/8] ^= (1 << ((bit_position + i) % 8));
            }

            // Analyse LOGIC
            int original_bit = (original_message[(bit_position + i)/8]>>((bit_position + i)%8))&1;
            if (original_bit == 0) flips += ones;
            else                   flips += (RHO - ones);
            freq[original_bit][ones]++;
        }
    }

    void MessageDecoder(unsigned char *plainText, unsigned char** message, unsigned char *original_message) {
        memset(plainText, 0, (K * L)/8);
        for(int i = 0; i < L; i++) {
            MessageDecoderP(plainText, K * i, message[i], original_message);
        }
    }
	
	void alloc2DCharArray(unsigned char** &arr) {
		arr = new unsigned char*[L];
		for(int i = 0; i < L; i++) {
            arr[i] = new unsigned char[N/8 + 1];
        }
	}
	
	void dealloc2DCharArray(unsigned char** &arr) {
		for(int i = 0; i < L; i++) 
			delete[] arr[i];
		delete[] arr;
	}
	
	void allocArray() {
		HP_bytes = new unsigned char[N/8 + 1];
		
		for(int i = 0; i < 2; i++) {
            freq[i] = new double[RHO + 1];
            for(int j = 0; j <= RHO; j++)
                freq[i][j] = 0;
        }
        
        int message_len = K * L;
        message = new unsigned char[(message_len + 7) / 8];
        recovered_message = new unsigned char[(message_len + 7) / 8];
        
        alloc2DCharArray(c2);
        alloc2DCharArray(encoded_message);
        alloc2DCharArray(decoder_inp);
	}
	
	void deallocArray() {
		delete[] HP_bytes;
		
		delete[] freq[0];
        delete[] freq[1];
        
		delete[] message; 
        delete[] recovered_message;
       
        dealloc2DCharArray(c2);
        dealloc2DCharArray(encoded_message);
        dealloc2DCharArray(decoder_inp);
	}
	
    public:

    MersennePoly(int N, int K, int H, int L, int rand_seed): N(N), K(K), H(H), L(L), rand_seed(rand_seed),
		key_generator(N, H, L) {
    // N, K, H, L, seed
	
        SetSeed(ZZ(rand_seed));
        mersenne_twister.seed(rand_seed);

        ZZ p = power(ZZ(2), N) - 1;
        ZZ_p::init(p);

        ZZ_pX P;
        SetCoeff(P, L, 1);
        SetCoeff(P, 0, -1);

        ZZ_pE::init(P);
        
        RHO = N/K;
        allocArray();
    }
    ~MersennePoly() {
		deallocArray();
	}

    Pair<Pair<ZZ_pE, ZZ_pE>, ZZ_pE> GenerateKeys() {
        ZZ_pE f = GenerateHSparse();
        ZZ_pE g = GenerateHSparse();
        ZZ_pE r = GenerateRandom();
        ZZ_pE t = f*r + g;
        return cons(cons(r, t), f);
    }
    
    Pair<Pair<ZZ_pE, ZZ_pE>, ZZ_pE> GenerateKeys_SHAKE() {
		
        ZZ_pE f = GenerateHSparse_SHAKE();
        ZZ_pE g = GenerateHSparse_SHAKE();
        ZZ_pE r = GenerateRandom_SHAKE();
        ZZ_pE t = f*r + g;
        return cons(cons(r, t), f);
    }

    void Encrypt(ZZ_pE& c1, unsigned char** c2, Pair<ZZ_pE, ZZ_pE> public_key, unsigned char *message) {
        ZZ_pE r = public_key.a, t = public_key.b;
        ZZ_pE a = GenerateHSparse();
        ZZ_pE b1 = GenerateHSparse(); 
        ZZ_pE b2 = GenerateHSparse();
        c1 = r*a + b1;
        MessageEncoder(encoded_message, message);
        XOR(c2, (t*a + b2), encoded_message);
        // (c1, c2) = (ra + b1, (ta + b2) xor encoded_message) 
    }

    void Decrypt(unsigned char *plainText, ZZ_pE private_key, ZZ_pE c1, unsigned char** c2, unsigned char *original_message) {
        
        XOR(decoder_inp, private_key * c1, c2);
        MessageDecoder(plainText, decoder_inp, original_message);
        // run message_decoder(fc1 xor c2)
    }

   vector<double> runs_final(int nruns, const char* fileName) {
        
        flips = 0;

        int pass = 0;
        auto t1 = chrono::high_resolution_clock::now();
        
        for(int big = 0; big < nruns; big++) {
            Pair<Pair<ZZ_pE, ZZ_pE>, ZZ_pE> keys = GenerateKeys_SHAKE();
            Pair<ZZ_pE, ZZ_pE> public_key = keys.a;
            ZZ_pE private_key = keys.b;

            
            for(int i = 0; i < (K * L) / 8; i++) {
                message[i] = RandomBnd(256);
            }
			
            Encrypt(c1, c2, public_key, message); /// c1, c2 form cipherText 
            Decrypt(recovered_message, private_key, c1, c2, message);

            bool ok = 1;
            for(int i = 0; i < (K * L)/8; i++) {
                if(message[i] != recovered_message[i]) {
                    ok = 0;
                }
            }
            if(ok)  pass++; 
        }
        auto t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> total_ms = t2 - t1;

        ofstream myFile;
        myFile.open(fileName, ofstream::out | ofstream::app);
        myFile << N << ", " << K << ", " << H << ", " << L << ", " << nruns << ", " << rand_seed << endl;

        double total = 0;
        for(int i = 0; i <= 1; i++)
            for(int j = 0; j <= RHO; j++)
                total += freq[i][j];

        double mean0 = 0, mean1 = 0;
        double var0 = 0, var1 = 0;

        for(int i = 0; i <= 1; i++) {
            double mean = 0, variance = 0;
            double local_total = 0;
            for(int j = 0; j <= RHO; j++) {
                mean += freq[i][j] * j;
                local_total += freq[i][j];
            }
            mean /= local_total;
            for(int j = 0; j <= RHO; j++) {
                variance += freq[i][j] * (j - mean)*(j - mean);
            }
            variance /= double(local_total - 1);

            for(int j = 0; j < RHO; j++) {
                freq[i][j] /= double(total);
                myFile << freq[i][j] << ", ";
                // freq[0][...] row
                // freq[1][...] row
            }
            freq[i][RHO] /= double(total);
            myFile << freq[i][RHO] << endl;
            myFile << mean << ", " << variance << endl;
            // mean and variance for original bit = 0
            // mean and variance for original bit = 1
            if(i == 0) {
                mean0 = mean;
                var0 = variance;
            } else {
                mean1 = mean;
                var1 = variance;
            }
        }
        myFile.close();
		
		double accuracy = double(pass) / double(nruns);
        double avg_ms = total_ms.count() / double(nruns);
        double erfc1 = (double(RHO / 2) - mean0) / sqrt(2 * var0);
        double erfc2 = (mean1 - double(RHO / 2)) / sqrt(2 * var1);
        double avg_flips = flips / double(nruns);
        return {accuracy, avg_ms, erfc1, erfc2, avg_flips};
    }
};
