for fidelity distance use

! F=tr { p**0.5 * sigma * p**0.5 }**0.5

! 1 diag p: p=A*D*A**H
! p**0.5 = A* D**0.5 * A**H

! 2. change basis of sigma to match new p,
! A * sigma * A**H = B

! subs B in
! F = tr{ A * D**0.5 * A**H) * B * (A * D**0.5 * A**H) }**0.5

! now F = tr{ A * D**0.5 * A**H * (A * sigma * A**H) * A * D**0.5 * A**H }**0.5

! F = tr{ A * D**0.5 * sigma_D * D**0.5 * A**H }**0.5

! F = tr{ A * C * A**H } **0.5

