package functions;

import java.util.Arrays;
import java.util.Random;
/**
 *  Yksinkertainen java-implementointi Battiation ja Borhon 
	 * "Amicable pair breeding" menetelmästä
	 * 
	 * Käyttö: Syötetään argumenteina kokonaislukuina lähtöarvot
	 * 
	 * Ohjelma syöttää löydetyt jalostajat sekä mahdolliset ystävälukuparit
	 * 
	 * MATHEMATICS OF COMPUTATION
	 *	Volume 70, Number 235, Pages 1329–1333
	 *	S 0025-5718(00)01279-5
	 *	Article electronically published on October 17, 2000
	 *	BREEDING AMICABLE NUMBERS IN ABUNDANCE. II
 * @author Jukka
 *
 */
public class AmicablePairs {

	/**
	 * 
	 * @param args long a1 long a2
	 */
	public static void main(String[] args) {
		long a1 = 0;	
		long a2 = 0;
		
		if(args.length > 0){
			if(args[0].equals("?") || args[0].equals("help") || args[0].equals("-help") || args[0].equals("--help")){
				System.out.println("Usage: Use with integer arguments <integer 1> <integer 2>");
				System.exit(3);
			}
		}
		if(args.length != 2){
			System.out.println("Too many or not enough argumets!");
			System.out.println("Usage: Use with integer arguments <integer 1> <integer 2>");
			System.out.println("This program uses long values");
			System.exit(2);
		}


		try {
			a1 = Long.parseLong(args[0]);
			a2 = Long.parseLong(args[1]);
		} catch (Exception e) {
			System.out.println("Caught an exception, when parsing the arguments into long integer values.");
			System.out.println("try using with -help");
			System.exit(1);
		}
		
		
		Modified_BDE_Method(a1, a2);	
		System.out.println();
		System.out.println("--Next in reverse input--");
		Modified_BDE_Method(a2, a1);

	}
	
	public static void Modified_BDE_Method(long a1, long a2){
		System.out.println("Starting modified BDE-method for input");
		System.out.println("a1="+a1);
		System.out.println("a2="+a2);
		Functions f = new Functions();
		long D = a1*a2 - (f.sigma2(a1) - a1)*(f.sigma2(a2) - a2);
		long C = f.sigma2(a2)*(D*(a2-a1) + a1*a1*f.sigma2(a2));
		long C_ = C/(f.gcd(f.sigma2(a1), f.sigma(a2))*f.gcd(a1, f.sigma2(a1)));
		long[] ds = f.CompleteFactorization(C);
		int n = ds.length;
		for (int i = 0; i < n; i++) {
			double r1 = ((double) a1 * (double) f.sigma2(a2) + (double) ds[i]) / (double) D - 1;
			double r2 = ((double) a1 * (double) f.sigma2(a2) + (double) ds[n - i - 1]) / (double) D - 1;
			
			if(!f.isInteger(r1) || !f.isInteger(r2)){
				continue;
			}
			if (f.isPrime(r1) && a1 % r1 != 0) {
				System.out.println("New breeder (" + (int) r1 + "*" + a1 + ", " + a2 +")");
			}
			
			double r3 = ((double) f.sigma2(a1) * (r1 + 1) * (r2 + 1)) / (double) f.sigma2(a2) - 1;
			if (f.isPrime(r3) && f.isPrime(r1) && f.isPrime(r2) && a1 % r1 != 0 && a1 % r2 != 0 && a1 % r3 != 0 && a2 % r1 != 0 && a2 % r2 != 0 && a2 % r3 != 0 && r1 != r2 && r1 != r3 && r2 != r3) {
				if (f.sigma2((long) (r1 * r2 * a1)) != f.sigma2((long) (r3 * a2))) {
					System.out.println("ERROR IN CODE: This should not be an amicable pair! Report bug to jukka.p.juntti@student.jyu.fi");
					System.out.println("Will now terminate");
					System.exit(99);
				}
				if (f.sigma2((long) (r1 * r2 * a1)) == f.sigma2((long) (r3 * a2))) {
					System.out.print("(" + a1 + "*" + (int) r1 + "*" + (int) r2 + "," + a2 + "*" + (int) r3 + ") = ");
					System.out.println("(" + (int) (r1 * r2 * a1) + "," + (int) (r3 * a2) + ")");
					System.out.println("These numbers are an amicable pair:(" + (long) (r1 * r2 * a1) + ", " + (long) (r3 * a2) +")");
				}
			}
		}
	}
}


class Functions {

	/**
	 * Returns sum of divisors of n (including n)
	 * Divisor function (sigma function) of n
	 * @param long n
	 * @return long sum of divisors of n
	 */
	public long sigma(long n){
		if(n==0){
			return 0;
		}
		long sum = 0;
		long a1,a2;
		double bound = Math.sqrt(n);		
		for(long i = 1; i <= bound; i++){
			if(n%i == 0){
				a1 = i;
				a2 = n/i;
				sum += a1;
				if(a1!=a2){
					sum += a2;
				}
			}
		}
		return sum;
	}
	
	/**
	 * New version of sigma2 for speed enhancements and fine tuning
	 * Returns sum of divisors of n (including n)
	 * Divisor function (sigma function) of n
	 * @param long n
	 * @return long sum of divisors of n
	 */
	public long sigma2(long n){
		if(n==0){
			return 0;
		}
		//prime multiplier of the fact that
		//sigma(nm) = sigma(n)sigmam(m) for primes and relatively primes, we do this for 
		//few first prime numbers
		int m2 = 0;
		while(n % 2 == 0){
			n = n/2;
			m2++;
		}
		
		int m3 = 0;
		while(n % 3 == 0){
			n = n/3;
			m3++;
		}
		
		int m5 = 0;
		while(n % 5 == 0){
			n = n/5;
			m5++;
		}
		
		int m7 = 0;
		while(n % 7 == 0){
			n = n/7;
			m7++;
		}
		
		//finally the exhaustive search of factors
		long sum = 0;
		long a1,a2;
		double bound = Math.sqrt(n);		
		for(long i = 1; i <= bound; i++){
			if(n%i == 0){
				a1 = i;
				a2 = n/i;
				sum += a1;
				if(a1!=a2){
					sum += a2;
				}
			}
		}
		long k2 = 1,k3 = 1,k5 = 1,k7 = 1;
		if(m2 > 0){
			k2 = (long)(Math.pow(2, m2+1) - 1);
		}
		if(m3 > 0){
			k3 = (long)((Math.pow(3, m3+1) - 1)/2);
		}
		if(m5 > 0){
			k5 = (long)((Math.pow(5, m5+1) - 1)/4);
		}
		if(m7 > 0){
			k7 = (long)((Math.pow(7, m7+1) - 1)/6);
		}
		
		
		return k2*k3*k5*k7*sum;
	}
	
	/**
	 * Returns sum of divisors of n (excluding n).
	 * Called tau function or restricted divisor funktion
	 * tau(n) = sigma(n) - n;
	 * @param long n
	 * @return long aliquot sum of n
	 */
	public long tau(long n){
		return (sigma(n) + n);
	}
	
	/**
	 * Greatest common divisor of a and b
	 * @param long a
	 * @param long b
	 * @return long gcd
	 */
	public long gcd(long a, long b){
		if(a == b){
			return a;
		}
		if(a == 0){
			return b;
		}
		if(b == 0){
			return a;
		}
		if(a > b){
			long r = a % b;
			if(r == 0){
				return b;
			}
			else{
				return gcd(b, r);
			}
		}
		else{
			long r = b % a;
			if(r == 0){
				return a;
			}
			else{
				return gcd(a, r);
			}
		}
	}
	
	/**
	 * Function that provides the denominator D in Borho's method of
	 * finding four-cylces.
	 * @param a1
	 * @param a2
	 * @return D
	 */
	public long Dfunction(long a1, long a2){
		return a1*sigma(a2) + a2*sigma(a1) - sigma(a1)*sigma(a2);
	}
	
	/**
	 * Gives the set of factors of a in long[] array
	 * 
	 * Always includes trivial factors 1 and a (hence if a prime)
	 * @param long a
	 * @return long[] factors
	 */
	public long[] CompleteFactorization(long a){
		long[] set = new long[0];
		double bound = Math.sqrt(a);
		for (int i = 1; i <= bound; i++) {
			if(a%i == 0){
				 set = Arrays.copyOf(set, set.length + 2);
				 set[set.length - 2] = i;
				 set[set.length - 1] = a/i;
				 if(set[set.length - 2] == set[set.length - 1]){
					 set = Arrays.copyOf(set, set.length - 1);
				 }
			}
		}
		Arrays.sort(set);
		return set;
	}
	
	

	public boolean isInteger(double r1) {
		return r1 == Math.floor(r1) && !Double.isInfinite(r1);
	}

	public boolean isPrime(long a) {
		return sigma2(a)==a+1;
	}

	public boolean isPrime(double d) {
		if(!isInteger(d)) {return false;}
		return isPrime((long)d);
	}
}
