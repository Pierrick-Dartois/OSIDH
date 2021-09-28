This is the source code of our implementation of OSIDH due to Colò and Kohel 
(2019) and our attack in Python/Sagemath.

(C) Anonymous authors, 2021.


Content of the main directory
=====================================================================

The source code is made of the following files:

* Group_basis.py contains generic functions to deal with finite abelian
groups of smooth cardinality, especially to compute discrete logarithms 
and group basis. The implementation is a simplified version of Sutherland's 
algorithms (2010). Most of thes functions are not used in our attack 
because the ideal class groups we work with are very simple (cyclic or 
almost cyclic).

* Relations_lattice.py contains code to computes the relations lattice of 
the prime ideals lying above the t first splitting primes q_1,...,q_t in 
the ideal class group of the order O_n:=Z+l^n O_K, where K:=Q(i), l=2, n=256 
and t=74 (parameters of Colò and Kohel (2019), Section 6). Of course, the 
values of the parameters can be changed. The code uses functions of 
Group_basis.py and saves the result in an external file (in the directory
Data_files).

* Find_SVP.py contains code to find short vectors in the relations lattice 
using fpylll implementation of the BKZ algorithm. It uses the output file of 
Relations_lattice.py and saves the result in another external file.

* OSIDH_protocol.py contains our implementation of OSIDH. It is made of 
three classes with appropriate methods:
..1. OSIDH, to instanciate all the necessary data to run a protocol (e.g. 
prime ideals, modular polynomials...);
..2. Chain, to instanciate descending l-isogeny chains;
..3. Chain_hor, to instanciate horizontal chains used in the protocol (given
by the action of the prime ideals q_j).
The procedure OSIDH_exe executes the protocol when OSIDH is instanciated and
when a public descending l-isogeny chain is given. OSIDH_simple_exe executes
the naive broken Diffie-Hellman protocol proposed by Colò and Kohel (2019)
in Section 5.1.

* OSIDH_attack.py contains our implementation of our attack in three steps:
..1. Recovering the whole descending l-isogeny chains of Alice and Bob;
..2. Recovering Alice's secret ideal class exponents with the knowledge
of Alice's chain and the public chain;
..3. Recovering the secret shared chain by acting on Bob's chain with 
Alice's exponents. 
The procedure OSIDH_attack_exe executes the protocol and the attack when 
OSIDH is instanciated and when a public descending l-isogeny chain is given. 
The procedure OSIDH_attack_stats runs the protocol and the attack several 
times to provide execution time statistics.

* OSIDH_attack_tests.py tests the attack (OSIDH_attack_exe) for different
sets of toy parameters. 

* OSIDH_attack_stats.py uses the function OSIDH_attack_stats from 
OSIDH_attack.py to provide execution time statistics. The results are saved
in an external .csv file in the directory Data_files.

The main directory also contains two subdirectories:

* Modular_polynomials containing all modular polynomials of prime indices up 
to 101. Those were extracted from Sutherland's database https://math.mit.edu/
~drew/ClassicalModPolys.html.

* Data_files containing execution results. It already contains execution files:
..1. Relations_lattice_basis.txt is an output from Relations_lattice.py;
..2. Lattice_BKZ_4.txt is an output from Find_SVP.py;
..3. Attack_stats.csv is an output from OSIDH_attack_stats.py.
Note that the lattice basis matrix in Relations_lattice_basis.txt is represented 
in columns while the reduced BKZ matrix in Lattice_BKZ_4.txt is represented in 
rows.


How to run the code?
=====================================================================

You should have Sage installed on your computer. You can as well use the
online Sage environment on https://cocalc.com.

First, open a Sage terminal and cd the main directory with the source code
(OSIDH_submit if you did not rename it). Then, you can load the .py file you
want to run (instruction load("filename")).

If you want to run the complete attack with toy parameters
---------------------------------------------------------------------

To test the attack, run the file OSIDH_attack_tests.py.

To measure the execution time on multiple executions, run the file 
OSIDH_attack_stats.py.


If you want to test lattice reduction with realistic parameters
---------------------------------------------------------------------

Run Relations_lattice.py and then Find_SVP.py


Warnings on the choice of parameters
=====================================================================

The OSIDH instanciation is limited by the available modular polynomials
in the Modular_polynomials directory. For that reason, the number of primes
t is limited to 10 (or slightly above). 

To make sure the attack is possible and runs in reasonable time, one should
make sure (2r+1)^t is close to l^n so that the key space covers the ideal
class group without many redundancies. 
