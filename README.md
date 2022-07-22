
# Introduction #

This project is an C++ (but in C-sytle) implementation of the lattice enumeration with discrete pruning for solving SVP, 
according to the paper *Lattice Enumeration with Discrete Pruning: Improvement, Cost Estimation and Optimization*
<!-- The full version is available at ? -->

Author: 
Luan Luan

Henan Key Laboratory of Network Cryptography Technology, 
Zhengzhou City, Henan Prov, China.



# Requirements #

This C++ program is developed and tested under Linux (Ubuntu >= 18.04 is recomended).
Please ensure that you have already installed the following libraries,
which are required to support. 
The experimental result in our paper is tested under these given version:

* gcc-7.5.0 (https://gcc.gnu.org/)
* gmp-6.1.0 (https://gmplib.org/)
* mpfr-4.0.2 (https://www.mpfr.org/)
* NTL-11.3.2 (https://www.shoup.net/ntl/)
* fplll-5.3.1 (https://github.com/fplll/fplll)

Noting that you may be able to use NTL __older than 9.4__ to submit your records to SVP Challenge as mentioned in https://www.latticechallenge.org/svp-challenge/. 
You can also modify `main.cpp` to input a certain basis from external file (But I'm not sure what will happen if the lattice is not a `random` one).


# Compile and Usage #

* Go into the `/src/` folder in shell and run `make` to compile.

* Go back to the root directory and run 

```
/.DPENUM --dim <dimension of SVP challenge> \
         --seed <seed for gnerate lattice basis> \
         --M <parameter for cell enumeration> \
         --bkz <blocksize of BKZ for processing basis>\
         --tours <number of tours for reprocessing> \
         --approx <approximate factor for SVP target length = approx*GH(L)> \
         --simtype <simulator type (more detail refers to depenum.cpp)> \
         --insertion <use inserting techniques in FK algorithm> 
```

All the parameters have default value in `main.cpp`.

For example, run `./DPENUM --dim 80 --M 50000 --bkz 32` to run a random SVP challenge instance of dimension 80 with DP enumeration parameters `M=50000` and `beta=32`. Now only `dim<200` is supported.

If you input `--simtype 0/1/2`, the program will ignore the `M` and `bkz` parameter and find the optimal parameter for this instance automatically. (See `depenum.cpp` to get the detailed meaning of this parameter)

The `--insertion 1` means to use the technique of FK algorithm, which use some short lattice vectors to help reduce the basis. This function maybe can heuristically accelerate the algorithm, but more details is under testing now.


# Additional Notes #

We thank github user *Harry Liang* for his generous help in modifying the erfz function for our special needs. (https://github.com/lhprojects/erfz)

This project is only a single-thread implementation without any parallel optimization, therefore it doesn't achieve an ideal efficiency for now. 
And it only supports lattice dimension `n<200` because of the issue about precision setting. 
In the future we will working in these aspects:

* The massive parallelization of DP enumeration.

* Modify its precision setting to extend the experiment to higher dimension.

* Add a GS sequence simulator for the case that `\bata \sim O(n)`. Shi Bai gives a good rectified BKZ simulator (https://github.com/BKZsimulator/probabilistic_simulator) in sage, and we will consider to translate it into C++.

* Get rid of NTL lib. Actually this program doen't need NTL, but we use some data structures and functions of NTL lib for some issues left over from history. They are just so convenient.

The file `default.json` in root directory is a part of fplll library which is used by BKZ 2.0 reduction. Please don't remove it.

If you have any questions/ideas or meet with bugs, please contact lunaluan9555@gmail.com. I would very like to help but I'm not sure if I have enough time.