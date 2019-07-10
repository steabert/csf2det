using Combinatorics

function expand_csf(stepvector::String, twoms::Integer)
  big_output = ""
  # count number of molecular orbitals n_mo
  n_mo = length(stepvector)
  n_somo = count(c->(c=='u' || c=='d'), stepvector)
  n_domo = count(c->c=='2', stepvector)

  # Generate Paldus tableau
  a = zeros(Int64, n_mo)
  b = zeros(Int64, n_mo)
  c = zeros(Int64, n_mo)
  ch = stepvector[1]
  if ch=='0'
      c[1]+=1
  elseif ch=='u'
      b[1]+=1
  elseif ch=='d'
      a[1]+=1
      b[1]-=1
      c[1]+=1
  elseif ch=='2'
      a[1]+=1
  else
      error("Invalid stepvector element: ", ch)
  end

  for i in 2:n_mo
    a[i] = a[i-1]
    b[i] = b[i-1]
    c[i] = c[i-1]
    ch = stepvector[i]
    if ch=='0'
        c[i]+=1
    elseif ch=='u'
        b[i]+=1
    elseif ch=='d'
        a[i]+=1
        b[i]-=1
        c[i]+=1
    elseif ch=='2'
        a[i]+=1
    else
        error("Invalid stepvector element: ", ch)
    end
  end

  # number of electrons is doubly occ + ud couples + excess alpha 
  n_e = 2 * a[n_mo] + b[n_mo]
  # total spin in units of 1/2 is number of excess alpha 
  spin = b[n_mo]
  # number of possible alpha spins 
  n_alpha = Int((n_somo + twoms)/2)

  # begin long loop over combinations n_alpha out of n_somo 
  #println("$n_e electrons in $n_mo orbitals")
  #println("output = phase * C^2 * SD")
  det = Array{Char}(undef,n_mo)

  for alpha_somos in combinations(1:n_somo, n_alpha)
    # initialize the incremental variables
    coeff = 1
    i_domo = i_somo = 1;
    i_alpha = i_beta = 0;

    for i in 1:n_mo
      ch=stepvector[i]
      if ch=='0'
        det[i] = '0';
      elseif ch=='u'
        if i_somo in alpha_somos
          det[i] = 'a'
          coeff *= (a[i] + b[i] - i_beta)
          i_alpha+=1
        else
          det[i] = 'b'
          coeff *= (a[i] + b[i] - i_alpha)
          i_beta+=1
        end
        coeff //= b[i];
        i_somo+=1
      elseif ch=='d'
        if i_somo in alpha_somos
          det[i] = 'a'
          coeff *= (i_beta - a[i] + 1)
          i_alpha+=1
          if iseven(b[i])
              coeff *= -1
          end
        else
          det[i] = 'b'
          coeff *= (i_alpha - a[i] + 1)
          i_beta+=1
          if isodd(b[i])
              coeff *= -1
          end
        end
        coeff //= (b[i] + 2)
        i_somo+=1
      elseif ch=='2'
        det[i] = '2'
        if isodd(b[i])
            coeff *= -1
        end
        i_alpha+=1
        i_beta+=1
      else
          error("Invalid stepvector elementi: ", ch)
      end
    end

    if coeff!=0 # if coefficient different from 0
      output = string(coeff)
      output *= " |"
      output *= join(det)
      output *= "|\n"
      #println(output)
      big_output *= output
    end
  end
  return big_output
end
