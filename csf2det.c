/* csf2det.c */

/***************************************************************************************************
 * program csf2det
 * 
 * expand a configuration state function (csf)
 * into a linear combination of determinants
 * using the GUGA approach.
 *
 * written by
 *  Steven Vancoillie
 *  Willem Van den Heuvel
 * (January 2009)
 *
 * References:
 * -----------
 * Isaiah Shavitt, "GUGA and its applications to direct CI calculations",
 * in "The Unitary Group for the Evaluation of Electronic Matrix Elements",
 * edited by J. Hinze, Lecture Notes in Chemistry 22, Springer-Verlag, Berlin, 1981,p.55 
 *
 * Libraries:
 * ----------
 *  Argtable (command line parser)
 *    http://argtable.sourceforge.net/
 *
 **************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <argtable2.h>

struct arg_lit *help, *verbose;
struct arg_str *stepvec;
struct arg_int *twoms;
struct arg_end *end;

int csf2det (const char * stepvector, int twoms);

int
main(int argc, char* argv[])
{
    /* command line parsing through argtable package */
    void* argtable[] = {
        help        = arg_lit0(NULL, "help", "display this help and exit"),
        verbose     = arg_lit0("v", "verbose", "print extra information"),
        stepvec     = arg_str1("s", "stepvec", "\"{0,u,d,2}\"", "string with spin coupling, e.g. \"2udu u0\""),
        twoms       = arg_int1("m", "twoms", "2*Ms", "Ms in units of one half"),
        end         = arg_end(20) };
    const char* progname = "csf2det";
    int rc = 0;
    int nerrors;

    /* verify the argtable[] entries were allocated sucessfully */
    if (arg_nullcheck(argtable) != 0)
        {
        /* NULL entries were detected, some allocations must have failed */
        printf("%s: insufficient memory\n",progname);
        rc=1;
        goto exit;
        }

    /* set default values */

    /* parse the command line flags, overriding default values */
    nerrors = arg_parse(argc,argv,argtable);

    /* special case: '--help' takes precedence over error reporting */
    if (help->count > 0)
        {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout,argtable,"\n");
        printf("\n");
        arg_print_glossary(stdout,argtable,"  %-40s %s\n");
        printf("\n");
        rc=0;
        goto exit;
        }

    /* special case: no command line options induces brief help */
    if (argc==1) {
        printf("Try '%s --help' for more information.\n",progname);
        rc=0;
        goto exit;
    }

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0)
        {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,progname);
        printf("Try '%s --help' for more information.\n",progname);
        rc=1;
        goto exit;
        }

    /* set global structures */

    /* call csf2det program */
    rc = csf2det (stepvec->sval[0], twoms->ival[0]);

exit:
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    return rc;
}

struct kCn {
  int n;
  int k;
  int * lex;
};

struct frac {
  int nom;
  int denom;
};

int combination (struct kCn * comb);
int gcd(int a, int b);
int simplify (int * frac);

int
csf2det (const char * stepvec, int twoms)
{
  int rc = 0;
  int i;

  /* count number of molecular orbitals n_mo */
  int n_mo = 0;
  for (i = 0; i < strlen (stepvec); i++ ) {
    if (stepvec[i] != ' ') ++n_mo;
  }

  /* csf array stores Shavitt numbers */
  int * csf = malloc (n_mo * sizeof (int));

  /* count number of singly and doubly occupied mo's */
  int iorb = 0;
  int n_somo = 0;
  int n_domo = 0;
  for (i = 0; i < strlen (stepvec); i++ ) {
    switch (stepvec[i])
    {
      case '0' : csf[iorb++] = 0; break;
      case 'u' : csf[iorb++] = 1; ++n_somo; break;
      case 'd' : csf[iorb++] = 2; ++n_somo; break;
      case '2' : csf[iorb++] = 3; ++n_domo; break;
      case ' ' : break;
      default :
                 printf ("input error: illegal character \'%c\' in stepvector\n", stepvec[i]);
                 printf ("             check the -s or --stepvec input string\n");
                 return 1;
    }
  }

  /* generate Paldus a, b, and c arrays */
  int * a = malloc (n_mo * sizeof (int));
  int * b = malloc (n_mo * sizeof (int));
  int * c = malloc (n_mo * sizeof (int));
  a[0] = 0;
  b[0] = 0;
  c[0] = 0;
  switch (csf[0])
  {
    case 0 : ++c[0]; break;
    case 1 : ++b[0]; break;
    case 2 : ++a[0]; --b[0]; ++c[0]; break;
    case 3 : ++a[0]; break;
  }
  for (i = 1; i < n_mo; i++) {
    a[i] = a[i-1];
    b[i] = b[i-1];
    c[i] = c[i-1];
    switch (csf[i])
    {
      case 0 : ++c[i]; break;
      case 1 : ++b[i]; break;
      case 2 : ++a[i]; --b[i]; ++c[i]; break;
      case 3 : ++a[i]; break;
    }
  }

  /* check the b array for negative values:
   * these should not occur as this measures incremental total spin */
  for (i = 0; i < n_mo; i++) {
    if (b[i] < 0) {
      printf ("input error: invalid ud ordering in stepvector\n");
      printf ("             check the -s or --stepvec input string\n");
      return 1;
    }
  }

  /* number of electrons is doubly occ + ud couples + excess alpha */
  int n_e = 2 * a[n_mo-1] + b[n_mo-1];
  /* total spin in units of 1/2 is number of excess alpha */
  int spin = b[n_mo-1];

  /* check absolut twoms input, should be smaller than total spin */
  if (abs(twoms) > spin) {
    printf ("input error: exceeded maximum Ms value of\n");
    printf ("             -/+ %i half integer units\n", spin);
    printf ("             check the -m or --twoms input value\n");
    return 1;
  }
  /* Ms = -S, -S+2, ..., S-2, S in half integer units,
   *   so Ms should be odd or even, same as S
   */
  char * parity[2];
  parity[0] = "EVEN";
  parity[1] = "ODD";
  if ((spin + twoms) % 2) {
    printf ("input error: Ms should be an %s number of half integers\n", parity[(spin % 2)]);
    printf ("             check the -m or --twoms input value\n");
    return 1;
  }

  /* number of possible alpha spins */
  int n_alpha = (n_somo + twoms) / 2;

  /* Generate combinations of n_alpha out of n_somo:
   *   these are all possible determinants for this csf.
   * Then determine the coefficient for each determinant.
   *
   * First allocate determinant array (holds info 0, a, b, or 2 for each mo)
   * and allocate combination initialized to lexicographically first subset.
   * Then calculate the coefficient of that determinant
   * and proceed as long as there are combinations.
   */

  /* array holding 2, a, b, 0  used for printing */
  char * det = malloc (n_mo * sizeof (char));

  /* keep track of which somo is alpha in loop over determinants */
  struct kCn comb_alpha;
  comb_alpha.n = n_somo;
  comb_alpha.k = n_alpha;
  comb_alpha.lex = malloc (n_alpha * sizeof (int));

  int * spin_eigv = malloc (n_somo * sizeof (int));

  /* begin long loop over combinations n_alpha out of n_somo */
  int phase, frac[2];
  int i_domo, i_somo;
  int i_alpha, i_beta;
  printf ("%i electrons in %i orbitals\n", n_e, n_mo);
  printf ("output = phase * C^2 * SD\n");

  /* initialize first combination */
  for (i = 0; i < n_alpha; ++i) comb_alpha.lex[i] = i;

  do
  {
    /* define the spin_eigv vector */
    for (i = 0; i < n_somo; ++i) spin_eigv[i] = -1;
    for (i = 0; i < n_alpha; ++i) spin_eigv[comb_alpha.lex[i]] = 1;

    /* initialize the incremental variables */
    phase = frac[0] = frac[1] = 1;
    i_domo = i_somo = 0;
    i_alpha = i_beta = 0;

    for (i = 0; i < n_mo; ++i)
    {
      switch (csf[i])
      {
        case 0 :
          det[i] = '0';
          break;
        case 1 :
          if (spin_eigv[i_somo++] == 1) {
            det[i] = 'a';
            frac[0] *= (a[i] + b[i] - i_beta);
            ++i_alpha;
          } else {
            det[i] = 'b';
            frac[0] *= (a[i] + b[i] - i_alpha);
            ++i_beta;
          }
          frac[1] *= b[i];
          break;
        case 2 :
          if (spin_eigv[i_somo++] == 1) {
            det[i] = 'a';
            frac[0] *= (i_beta - a[i] + 1);
            ++i_alpha;
            if (!(b[i] % 2)) phase *= -1;
          } else {
            det[i] = 'b';
            frac[0] *= (i_alpha - a[i] + 1);
            ++i_beta;
            if (b[i] % 2) phase *= -1;
          }
          frac[1] *= (b[i] + 2);
          break;
        case 3 :
          det[i] = '2';
          if (b[i] % 2) phase *= -1;
          ++i_alpha;
          ++i_beta;
          break;
      }
      simplify (frac);
    }
    if (frac[0]) { /* if coefficient different from 0 */
      if (phase > 0) 
        printf (" %3c", '+');
      else
        printf (" %3c", '-');
      printf (" %3i/%-8i", frac[0], frac[1]);
      printf (" |");
      for (i = 0; i < n_mo; ++i) printf (" %c", det[i]);
      printf (" |");
      printf ("\n");
    }
  } while (combination (&comb_alpha));

  free (spin_eigv); free (comb_alpha.lex); free (det);
  free (a); free (b); free (c);
  free (csf);

  return rc;
}

int
gcd(int a, int b) 
{ 
     return ( b == 0 ? a : gcd(b, a % b) ); 
}

int
simplify (int * frac)
{
  if (frac[0] == 0) return 0;
  int div = gcd (frac[0], frac[1]);
  frac[0] /= div;
  frac[1] /= div;
  return div;
}

int
combination (struct kCn * comb) {
  int i;
  /* search backward to first possible increase */
  int ptr = comb->k - 1;
  while (ptr >= 0 && comb->lex[ptr] == comb->n - comb->k + ptr) --ptr;
  if (ptr >= 0) {
    ++comb->lex[ptr];
    for (i = 1; i < comb->k - ptr; ++i) comb->lex[ptr + i] = comb->lex[ptr] + i;
    return 1;
  } else {
    return 0;
  }
}
