use std::fmt;
use structopt::StructOpt;

#[derive(StructOpt)]
struct Cli {
    stepvector: String,
    twoms: i32,
}

fn main() {
    let args = Cli::from_args();
    let result = csf2det(&args.stepvector, args.twoms);
    match result {
        Ok(msg) => println!("{}", msg),
        Err(error) => eprintln!("Error: {}", error),
    }
}

//
// Math stuff
//

// Simplify fractions

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

struct Fraction {
    nom: i32,
    denom: i32,
}

impl Fraction {
    fn simplify(&self) -> Fraction {
        if self.nom == 0 {
            Fraction { nom: 0, denom: 0 }
        } else {
            let div = gcd(self.nom, self.denom);
            Fraction {
                nom: self.nom / div,
                denom: self.denom / div,
            }
        }
    }
}

impl fmt::Display for Fraction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.nom == 0 || self.denom == 1 {
            write!(f, "{}", self.nom)
        } else {
            write!(f, "{}/{}", self.nom, self.denom)
        }
    }
}

// Combinatorics

struct Combination {
    n: usize,        // size of the set
    k: usize,        // size of the subset of distinct elements
    lex: Vec<usize>, // lexicographic representation kCn
}

fn next_combination(c: &mut Combination) -> bool {
    if c.k == 0 {
        return false;
    }
    let mut ptr = c.k - 1;
    while c.lex[ptr] == c.n - c.k + ptr {
        if ptr == 0 {
            return false;
        }
        ptr -= 1;
    }
    c.lex[ptr] += 1;
    for i in 1..(c.k - ptr) {
        c.lex[ptr + i] = c.lex[ptr] + i;
    }
    return true;
}

//
// Expand CSF to determinants
//

fn csf2det(stepvec: &str, twoms: i32) -> Result<String, String> {
    let spin_projection = Fraction {
        nom: twoms,
        denom: 2,
    };
    println!("CSF: {}, Ms: {}", stepvec, spin_projection.simplify());

    let n_mo = stepvec.len();
    let mut n_somo: usize = 0;

    // construct Pldus tableau
    let mut a: Vec<i32> = Vec::new();
    let mut b: Vec<i32> = Vec::new();
    let mut c: Vec<i32> = Vec::new();
    for (i, step) in stepvec.chars().enumerate() {
        if i == 0 {
            a.push(0);
            b.push(0);
            c.push(0);
        } else {
            a.push(a[i - 1]);
            b.push(b[i - 1]);
            c.push(c[i - 1]);
        }

        if step == '0' {
            c[i] += 1;
        } else if step == 'u' {
            n_somo += 1;
            b[i] += 1;
        } else if step == 'd' {
            n_somo += 1;
            a[i] += 1;
            b[i] -= 1;
            c[i] += 1;
        } else if step == '2' {
            a[i] += 1;
        } else {
            return Err(format!("invalid step: {}", step));
        }
    }

    // number of electrons
    let n_e = 2 * a[n_mo - 1] + b[n_mo - 1];
    // total spin in units of 1/2
    let spin = b[n_mo - 1];
    if twoms.abs() > spin || (spin + twoms) % 2 != 0 {
        return Err(format!(
            "half-integer Ms {} does not exist for spin {}",
            twoms, spin
        ));
    }

    println!("Determinants for {} electrons in {} orbitals:", n_e, n_mo);

    // number of possible alpha spin
    let n_alpha = ((n_somo as i32 + twoms) / 2) as usize;
    let mut comb_alpha = Combination {
        n: n_somo,
        k: n_alpha,
        lex: (0..n_alpha).collect(),
    };

    loop {
        // Set up the phase and fractional coefficient
        let mut phase = 1;
        let mut frac = Fraction { nom: 1, denom: 1 };

        // Orbital reference indices
        let mut i_somo: usize = 0;

        // alpha/beta parity
        let mut i_alpha = 0;
        let mut i_beta = 0;

        // set up spin eigenvalue vector
        let mut spin_eigv: Vec<i8> = vec![-1; n_somo];
        for i in 0..n_alpha {
            spin_eigv[comb_alpha.lex[i]] = 1;
        }

        // Construct the determinant.
        let mut det: String = String::new();
        for (i, step) in stepvec.chars().enumerate() {
            if step == '0' {
                det.push('0');
            } else if step == 'u' {
                if spin_eigv[i_somo] == 1 {
                    det.push('a');
                    frac.nom *= a[i] + b[i] - i_beta;
                    i_alpha += 1;
                } else {
                    det.push('b');
                    frac.nom *= a[i] + b[i] - i_alpha;
                    i_beta += 1;
                }
                frac.denom *= b[i];
                i_somo += 1;
            } else if step == 'd' {
                if spin_eigv[i_somo] == 1 {
                    det.push('a');
                    frac.nom *= i_beta - a[i] + 1;
                    i_alpha += 1;
                    if (b[i] % 2) == 0 {
                        phase *= -1;
                    }
                } else {
                    det.push('b');
                    frac.nom *= i_alpha - a[i] + 1;
                    i_beta += 1;
                    if (b[i] % 2) != 0 {
                        phase *= -1;
                    }
                }
                frac.denom *= b[i] + 2;
                i_somo += 1;
            } else if step == '2' {
                det.push('2');
                if (b[i] % 2) != 0 {
                    phase *= -1;
                }
                i_alpha += 1;
                i_beta += 1;
            }
            frac = frac.simplify();
        }

        if frac.nom != 0 {
            let sign = if phase > 0 { '+' } else { '-' };
            println!("{} {} |{}|", sign, frac, det);
        }

        // Continue with the next combination, or quit
        if !next_combination(&mut comb_alpha) {
            break;
        }
    }

    // return code
    Ok(String::from("done!"))
}
