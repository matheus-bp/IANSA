function brent(f, a, b; tol=1e-8, maxiter=100)
    #=
    Brent's method for finding a root of f(x) = 0 in the interval [a, b]
    
    Args:
        f: Function whose root we want to find (must be continuous)
        a, b: Initial bracket [a, b] where f(a) and f(b) have opposite signs
        tol: Tolerance for convergence (default: 1e-8)
        maxiter: Maximum number of iterations (default: 100)
    
    Returns:
        The approximate root x where f(x) ≈ 0
    
    Throws:
        Error if root is not bracketed or if max iterations exceeded
    =#
    
    # Initial function evaluations
    fa = f(a)
    fb = f(b)
    
    # Check that root is bracketed (f(a) and f(b) must have opposite signs)
    fa * fb >= 0 && error("Root not bracketed: f(a) and f(b) must have opposite signs")
    
    # Ensure b is the best guess (closer to root)
    if abs(fa) < abs(fb)
        a, b = b, a  # Swap a and b
        fa, fb = fb, fa
    end
    
    # Initialize variables:
    c = a       # Previous value of b
    d = c       # Second previous value
    fc = fa     # f(c)
    mflag = true # Flag to control bisection
    
    for i = 1:maxiter
        # Check if we've converged
        if abs(b - a) < tol || abs(fb) < tol
            return b
        end
        
        # Try inverse quadratic interpolation if possible
        if fa != fc && fb != fc
            # Inverse quadratic interpolation formula:
            s = a*fb*fc/((fa-fb)*(fa-fc)) + 
                b*fa*fc/((fb-fa)*(fb-fc)) + 
                c*fa*fb/((fc-fa)*(fc-fb))
        else
            # Fall back to secant method
            s = b - fb*(b - a)/(fb - fa)
        end
        
        # Check if the new estimate is acceptable:
        # Condition 1: s must be between (3a + b)/4 and b
        cond1 = (s < (3a + b)/4) || (s > b)
        
        # Condition 2: If mflag is set and |s-b| ≥ |b-c|/2
        cond2 = mflag && (abs(s - b) >= abs(b - c)/2)
        
        # Condition 3: If mflag is clear and |s-b| ≥ |c-d|/2
        cond3 = !mflag && (abs(s - b) >= abs(c - d)/2)
        
        # Condition 4: If bisection was used last time and |b-c| < tol
        cond4 = mflag && (abs(b - c) < tol)
        
        # Condition 5: If interpolation was used last time and |c-d| < tol
        cond5 = !mflag && (abs(c - d) < tol)
        
        # If any conditions are met, use bisection instead
        if cond1 || cond2 || cond3 || cond4 || cond5
            s = (a + b)/2  # Bisection step
            mflag = true   # Indicate we used bisection
        else
            mflag = false  # Indicate we used interpolation
        end
        
        # Evaluate function at new point
        fs = f(s)
        d = c   # Store previous c
        c = b   # Update c to previous b
        fc = fb # Update f(c)
        
        # Update bracket:
        if fa * fs < 0
            b = s
            fb = fs
        else
            a = s
            fa = fs
        end
        
        # Ensure b remains the best estimate
        if abs(fa) < abs(fb)
            a, b = b, a
            fa, fb = fb, fa
        end
    end
    
    error("Brent's method failed to converge in $maxiter iterations")
end