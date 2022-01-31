# Crypto
 
## Mersenne crypto project

First, install NTL library (instruction available at https://libntl.org/)
 
Use driver_final.cpp to run the crypto system.

```
make driver_final
./driver_final {security_parameter} (one of {128, 192, 256})
```
This will produce 2 files "freq_output{security_parameter}.csv" and "res_output{security_parameter}.csv"


For plotting of the result, use make_gauss.py (requires pyplot). 

```
python3 make_gauss.py file_name (in .csv format) 
```


