// CSF to determinants

function gcd(a, b)
{
  return b === 0 ? a : gcd(b, a % b);
}

function simplify (frac)
{
  if (frac[0] === 0) return 0;
  const div = gcd(frac[0], frac[1]);
  frac[0] /= div;
  frac[1] /= div;
  return div;
}

function next_combination (comb) {
  /* search backward to first possible increase */
  let ptr = comb.k - 1;
  while (ptr >= 0 && comb.lex[ptr] === comb.n - comb.k + ptr) --ptr;
  if (ptr >= 0) {
    ++comb.lex[ptr];
    for (let i = 1; i < comb.k - ptr; ++i) comb.lex[ptr + i] = comb.lex[ptr] + i;
    return true;
  } else {
    return false;
  }
}

function csf2det(stepvector, twoms) {
  // count number of molecular orbitals n_mo
  const n_mo = stepvector.length;
  let n_somo = 0;
  let n_domo = 0;

  // convert stepvector string to array: 0,u,d,2 → 0,1,2,3
  const csf = new Int32Array(n_mo);
  for (let i = 0; i < n_mo; i += 1) {
    switch (stepvector[i]) {
    case '0' : csf[i] = 0; break;
    case 'u' : csf[i] = 1; ++n_somo; break;
    case 'd' : csf[i] = 2; ++n_somo; break;
    case '2' : csf[i] = 3; ++n_domo; break;
    default:
      throw new Error(`invalid stepvector string: ${stepvector[i]}`);
    }
  }

  // Generate Paldus tableau
  const a = new Int32Array(n_mo);
  const b = new Int32Array(n_mo);
  const c = new Int32Array(n_mo);
  a[0] = 0;
  b[0] = 0;
  c[0] = 0;
  switch (csf[0]) {
  case 0 : ++c[0]; break;
  case 1 : ++b[0]; break;
  case 2 : ++a[0]; --b[0]; ++c[0]; break;
  case 3 : ++a[0]; break;
  }
  for (let i = 1; i < n_mo; i++) {
    a[i] = a[i-1];
    b[i] = b[i-1];
    c[i] = c[i-1];
    switch (csf[i]) {
    case 0 : ++c[i]; break;
    case 1 : ++b[i]; break;
    case 2 : ++a[i]; --b[i]; ++c[i]; break;
    case 3 : ++a[i]; break;
    }
  }

  /* number of electrons is doubly occ + ud couples + excess alpha */
  const n_e = 2 * a[n_mo-1] + b[n_mo-1];
  /* total spin in units of 1/2 is number of excess alpha */
  const spin = b[n_mo-1];
  /* number of possible alpha spins */
  const n_alpha = Math.floor((n_somo + twoms) / 2);

  const comb_alpha = {
    'n': n_somo,
    'k': n_alpha,
    'lex': new Int32Array(n_alpha)
  };

  const spin_eigv = new Int32Array(n_somo);

  /* begin long loop over combinations n_alpha out of n_somo */
  let phase;
  const frac = new Int32Array(2);
  let i_domo, i_somo;
  let i_alpha, i_beta;
  console.log(`${n_e} electrons in ${n_mo} orbitals`);
  console.log('output = phase * C^2 * SD');
  const det = Array(n_mo);

  /* initialize first combination */
  for (let i = 0; i < n_alpha; ++i) comb_alpha.lex[i] = i;

  do
  {
    /* define the spin_eigv vector */
    for (let i = 0; i < n_somo; ++i) spin_eigv[i] = -1;
    for (let i = 0; i < n_alpha; ++i) spin_eigv[comb_alpha.lex[i]] = 1;

    /* initialize the incremental variables */
    phase = frac[0] = frac[1] = 1;
    i_domo = i_somo = 0;
    i_alpha = i_beta = 0;

    for (let i = 0; i < n_mo; ++i) {
      switch (csf[i]) {
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
    let output = '';
    if (frac[0]) { /* if coefficient different from 0 */
      if (phase > 0) {
        output += ' + ';
      } else {
        output += ' - ';
      }
      output += `${frac[0]}/${frac[1]}`;
      output += ' |';
      for (let i = 0; i < n_mo; ++i) output += det[i];
      output += ' |';
      console.log(output);
    }
  } while (next_combination(comb_alpha));
}

csf2det('22ud0ud02', 0);
