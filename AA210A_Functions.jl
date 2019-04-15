"""

-----------------------------------------

**AA210A Compressible Flow Formulae**

-----------------------------------------

Material from "AA210A Course Reader: Fundamentals of Compressible Flow"

by Prof. Brian J. Cantwell - https://web.stanford.edu/~cantwell/

|**Function Name**| **Function Description**                   |**Eqn. No.**|
|:----------------|:-------------------------------------------|:-----------|
|           gamma | Specific Heat Ratio                        | (2.98)     |
|             C_p | Heat Capacity at Constant Pressure, J/kg-K | (2.98)     |
|             C_v | Heat Capacity at Constant Volume, J/kg-K   | (2.98)     |
|        spdofsnd | Speed of Sound, m/s                        | (2.112)    |
|             PtP | Stagnation/Static Pressure Ratio           | (5.85)     |
|             TtT | Stagnation/Static Temperature Ratio        | (5.77)     |
|   **NormShock** | **Normal Shock Jump Relations**            | **⤓⤓⤓**  |
|     `Mach2Norm` | `Normal Shock Exit Mach`                   | `(9.60)`   |
|      `T2T1Norm` | `Normal Shock Static Temperature Ratio`    | `(9.73)`   |
|      `P2P1Norm` | `Normal Shock Static Pressure Ratio`       | `(9.56)`   |
|    `Pt2Pt1Norm` | `Normal Shock Stagnation Pressure Ratio`   | `(9.66)`   |
|      `r2r1Norm` | `Normal Shock Density Ratio`               | `(9.97)`   |
|        FindMach | Find Mach Number from Parameter Ratios     |            |
|    **OblShock** | **Oblique Shock Jump Relations**           | **⤓⤓⤓**  |
|       `deflang` | `Oblique Shock Flow Turn Angle θ (deg)`    | `(12.15)`  |
|    `deflangsol` | `Oblique Shock Angle Solver`               | `(12.15)`  |
|      `Mach2Obl` | `Oblique Shock Exit Mach`                  | `(12.7)`   |
|       `T2T1Obl` | `Oblique Shock Static Temperature Ratio`   | `(12.7)`   |
|       `P2P1Obl` | `Oblique Shock Static Pressure Ratio`      | `(12.7)`   |
|     `Pt2Pt1Obl` | `Oblique Shock Stagnation Pressure Ratio`  | `(12.8)`   |
|       `r2r1Obl` | `Oblique Shock Density Ratio`              | `(12.7)`   |
|      **aratio** | **A*/A = f(M)**                            |**(10.16)** |
|     `aratiosol` | `A*/A = f(M) Solver`                       | `(10.16)`  |
|           Fanno | Fanno Line Relations                       | (11.31)    |
|        Rayleigh | Rayleigh Line Relations                    | (11.47)    |
|           PMang | Prandtl-Meyer Expansion Angle ω (deg)      | (12.38)    |
|        PMangsol | Prandtl-Meyer Expansion Angle Solver       | (12.38)    |
|   ShockTubeP4P1 | Shock Tube Initial Pressure Ratio          | (13.64)    |
|     ShockTubeUp | Shock Tube Piston Velocity Up              | (13.59)    |
|      ChokedCond | Pressure Ratio for Choked Flow Pt/Pa       | (10.26)    |
|      NozzleMach | Fully Expanded Nozzle Exit Mach Number     | (10.25)    |
|         scaleht | Atmospheric Scale Height, m                | (2.119)    |
|         shockRe | Shock Reynolds Number (ρUδ)/μ              | (9.96)     |
|       CJDetMach | Chapman-Jouget Detonation Mach Number      | (11.61)    |
 
Type \"AA210AFnList()\" or \"?AA210A_Functions\" to see function list again

Type \"?<function_name>\" to see individual function documentation

"""
module AA210A_Functions

export spdofsnd,
PtP,TtT,
NormShock,Mach2Norm,T2T1Norm,P2P1Norm,Pt2Pt1Norm,r2r1Norm,
FindMach,
OblShock,deflang,T2T1Obl,P2P1Obl,Pt2Pt1Obl,r2r1Obl,Mach2Obl,deflangsol,
aratio,aratiosol,
Fanno,Rayleigh,
PMang,
PMangsol,
ShockTubeP4P1,ShockTubeMs,ShockTubeUp, #ShockTubeMrs,
gamma,C_p,C_v,
ChokedCond,NozzleMach,
scaleht,
shockRe,
CJDetMach,
AA210AFnList

begin #local scope to avoid naming variables used elsewhere

local fnnames = ["Function Name|Function Description|Eqn. Nos."]
#name format: fnname(args)|description|eqn number

#~~~~~~~~~~Constants~~~~~~~~~~

# println("Ru = 8314.46: Universal Gas Constant, J/kmol-K")
# push!(fnnames,"Ru = 8314.46|Universal Gas Constant, J/kmol-K|")
# Ru = 8314.46 #Universal Gas Constant, J/kmol K

# println("g0 = 9.80665: Std. Gravitational Acceleration, m/s²")
# push!(fnnames,"g0 = 9.80665|Std. Gravitational Accel., m/s²|")
# g0 = 9.80665 #Std. Gravitational Acceleration

# push!(fnnames,"||")

#~~~~~~~~~~Functions~~~~~~~~~~

push!(fnnames,"gamma|Specific Heat Ratio|(2.98)")
"""

### Syntax
    gamma(dof)
    gamma()

### Description
_Eqn. 2.98_

Calculates the specific heat ratio of an ideal gas with total number of molecular degrees of freedom 'dof'. Default corresponds to Air. Typical values for 'dof' are 5 for diatomic gases and 3 for monatomic gases (Section A.8).

### Arguments
'dof::Integer': Molecular degrees of freedom of gas (Default: 5)

### See Also
C_p, spdofsnd, scaleht

"""
function gamma(dof::Integer) #Specific Heat Ratio (eq 2.98)
  gamma = (dof+2)/dof
  return gamma
end
function gamma() #Alternative Method for generic case (Air, diatomic gases)
  gamma(5)
end

push!(fnnames,"C_p|Heat Capacity at Constant Pressure, J/kg-K|(2.98)")
"""

### Syntax
    C_p(dof, Mw)
    C_p()

### Description
_Eqn. 2.98_

Calculates the specific heat capacity Cp of an ideal gas with total number of molecular degrees of freedom 'dof' and molecular weight 'Mw'. Default corresponds to Air. Typical values for dof are 5 for diatomic gases and 3 for monatomic gases (Section A.8).

### Arguments
'dof::Integer': Molecular degrees of freedom of gas (Default: 5)

'Mw::Real': Molecular weight, amu or kg/mol (Default: 28.97)

### See Also
gamma, C_v, spdofsnd, scaleht

"""
function C_p(dof::Real, Mw::Real) #Heat Capacity at constant pressure (eq 2.98)
  Ru = 8314.46 #Universal Gas Constant, J/kmol K
  Cp = ((dof + 2)/2)*(Ru/Mw)
  return Cp
end
function C_p() #Alternative Method for generic case (Air)
  C_p(5,28.97)
end

push!(fnnames,"C_v|Heat Capacity at Constant Volume, J/kg-K|(2.98)")
"""

### Syntax
    C_v(dof, Mw)
    C_v()

### Description
_Eqn. 2.98_

Calculates the specific heat capacity Cv of an ideal gas with total number of molecular degrees of freedom 'dof' and molecular weight 'Mw'. Default corresponds to Air. Typical values for dof are 5 for diatomic gases and 3 for monatomic gases (Section A.8).

### Arguments
'dof::Integer': Molecular degrees of freedom of gas (Default: 5)

'Mw::Real': Molecular weight, amu or kg/mol (Default: 28.97)

### See Also
gamma, C_p, spdofsnd, scaleht

"""
function C_v(dof::Real,Mw::Real) #Heat Capacity at constant volume (eq 2.98)
  Ru = 8314.46 #Universal Gas Constant, J/kmol K
  Cv = (dof*Ru)/(2*Mw)
  return Cv
end
function C_v() #Alternative Method for generic case (Air)
  C_v(5,28.97)
end

push!(fnnames,"spdofsnd|Speed of Sound, m/s|(2.112)")
"""

### Syntax
    spdofsnd(gam,Mw,T)
    spdofsnd(T)

### Description
_Eqn. 2.112_

Calculates the speed of sound of an ideal gas with specific heat ratio 'gam' and molecular weight 'Mw' at the temperature 'T'. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'Mw::Real': Molecular Weight, amu or kg/mol (Default: 28.97)

'T::Real': Static Temperature, K

### See Also
scaleht

"""
function spdofsnd(gamma::Real,Mw::Real,Temp::Real) #Ideal Gas Speed of Sound
  Ru = 8314.46 #Universal Gas Constant, J/kmol K
  sos = sqrt((gamma*Ru*Temp)/Mw)
  return sos
end
function spdofsnd(Temp::Real) #Alternative Method for generic case (Air)
  spdofsnd(1.4,28.97,Temp)
end

push!(fnnames,"PtP|Stagnation/Static Pressure Ratio|(5.85)")
"""

### Syntax
    PtP(gam,M)
    PtP(M)

### Description
_Eqn. 5.85_

Calculates the ratio of stagnation pressure to static pressure in a flow of an ideal gas with specific heat ratio 'gam' at Mach Number 'M'. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Mach Number

### See Also
TtT, P2P1Norm, Pt2Pt1Norm, T2T1Norm

"""
function PtP(gamma::Real,Mach::Real) #eq 9.62
  expo = gamma/(gamma-1)
  PtP = TtT(gamma,Mach)^expo
  return PtP
end
function PtP(Mach::Real) #Alternative Method for generic case (Air/diatomic)
  PtP(1.4,Mach)
end

push!(fnnames,"TtT|Stagnation/Static Temperature Ratio|(5.77)")
"""

### Syntax
    TtT(gam,M)
    TtT(M)

### Description
_Eqn. 5.77_

Calculates the ratio of stagnation temperature to static temperature in a flow of an ideal gas with specific heat ratio 'gam' at Mach Number 'M'. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Mach Number

### See Also
PtP, T2T1Norm, P2P1Norm, Pt2Pt1Norm

"""
function TtT(gamma::Real,Mach::Real) #eq 5.77
  TtT = 1 + (gamma-1)/2*Mach^2
  return TtT
end
function TtT(Mach::Real) #Alternative Method for generic case (Air/diatomic)
  TtT(1.4,Mach)
end

push!(fnnames,"NormShock|Normal Shock Jump Relations|(9.56-97)")
"""

### Syntax
    NormShock(which,gam,M)
    NormShock(which,M)

### Description
Calculates the normal shock jump relation selected by 'which' for an ideal gas with specific heat ratio 'gam' at shock Mach Number 'M'. Default 'gam' corresponds to Air.

### Arguments
'which::String': Shock Jump Relation Selection 
- "M" -> Mach Number immediately after shock, 'M2' _Eqn. 9.60_
- "T" -> Static temperature ratio across shock, 'T2/T1' _Eqn. 9.73_
- "P" -> Static pressure ratio across shock, 'P2/P1' _Eqn. 9.56_
- "Pt" -> Stagnation pressure ratio across shock, 'Pt2/Pt1' _Eqn. 9.66_
- "r" -> Density ratio across shock, 'r2/r1' _Eqn. 9.97_
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Shock Mach Number (Mach Number immediately before shock)

### See Also
OblShock, Mach2Norm, T2T1Norm, P2P1Norm, Pt2Pt1Norm, r2r1Norm

"""
function NormShock(ver::String,gamma::Real,Mach::Real) #Consolidated Normal Shock Jump Relations
  if ver == "M"
    Mach2Norm(gamma,Mach)
  elseif ver == "T"
    T2T1Norm(gamma,Mach)
  elseif ver == "P"
    P2P1Norm(gamma,Mach)
  elseif ver == "Pt"
    Pt2Pt1Norm(gamma,Mach)
  elseif ver == "r"
    r2r1Norm(gamma,Mach)
  end
end
function NormShock(ver::String,Mach::Real) #Alternative Method for generic case (Air/diatomic)
  NormShock(ver,1.4,Mach)
end

push!(fnnames,"Mach2Norm|*Normal Shock Exit Mach|(9.60)")
"""

### Syntax
    Mach2Norm(gam,M)

### Description
_Eqn. 9.60_

Calculates the mach number immediately after a normal shock of Mach Number 'M' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

### See Also
NormShock, Mach2Obl, T2T1Norm, P2P1Norm, Pt2Pt1Norm, r2r1Norm

"""
function Mach2Norm(gamma::Real,Mach::Real) #eq 9.60
  numer = 1 + (gamma-1)/2*Mach^2
  denom = gamma*Mach^2 - (gamma-1)/2
  mach2 = sqrt(numer/denom)
  return mach2
end

push!(fnnames,"T2T1Norm|*Normal Shock Static Temperature Ratio|(9.73)")
"""

### Syntax
    T2T1Norm(gam,M)

### Description
_Eqn. 9.73_

Calculates the static temperature ratio across a normal shock of Mach Number 'M' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

### See Also
NormShock, T2T1Obl, Mach2Norm, P2P1Norm, Pt2Pt1Norm, r2r1Norm

"""
function T2T1Norm(gamma::Real,Mach::Real) #eq 9.73
  numer = (1 + (gamma-1)/2*Mach^2)*(gamma*Mach^2 - (gamma-1)/2)
  denom = ((gamma+1)/2)^2*Mach^2
  T2T1 = numer/denom
  return T2T1
end

push!(fnnames,"P2P1Norm|*Normal Shock Static Pressure Ratio|(9.56)")
"""

### Syntax
    P2P1Norm(gam,M)

### Description
_Eqn. 9.56_

Calculates the static pressure ratio across a normal shock of Mach Number 'M' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

### See Also
NormShock, P2P1Obl, Mach2Norm, T2T1Norm, Pt2Pt1Norm, r2r1Norm

"""
function P2P1Norm(gamma::Real,Mach::Real) #Static Pressure across normal shock (eq. 9.56)
  Pratio = ((2*gamma)/(gamma-1)*Mach^2 - 1)/((gamma+1)/(gamma-1))
  return Pratio
end

push!(fnnames,"Pt2Pt1Norm|*Normal Shock Stagnation Pressure Ratio|(9.66)")
"""

### Syntax
    Pt2Pt1Norm(gam,M)

### Description
_Eqn. 9.66_

Calculates the stagnation pressure ratio across a normal shock of Mach Number 'M' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

### See Also
NormShock, Pt2Pt1Obl, P2P1Norm, Mach2Norm, T2T1Norm, r2r1Norm

"""
function Pt2Pt1Norm(gamma::Real,Mach::Real) #eq 9.66
  numer1 = (gamma+1)/(gamma-1)
  denom1 = (2*gamma)/(gamma-1)*Mach^2 - 1
  part1 = (numer1/denom1)^(1/(gamma-1))
  numer2 = (gamma+1)/2*Mach^2
  denom2 = 1 + (gamma-1)/2*Mach^2
  part2 = (numer2/denom2)^(gamma/(gamma-1))
  Pt2Pt1 = part1*part2
  return Pt2Pt1
end

push!(fnnames,"r2r1Norm|*Normal Shock Density Ratio|(9.97)")
"""

### Syntax
    r2r1Norm(gam,M)

### Description
_Eqn. 9.97_

Calculates the density ratio across a normal shock of Mach Number 'M' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

### See Also
NormShock, r2r1Obl, Pt2Pt1Norm, P2P1Norm, Mach2Norm, T2T1Norm

"""
function r2r1Norm(gamma,Mach) #eq 9.97
  r2r1 = ((gamma+1)*Mach^2)/((gamma-1)*Mach^2 + 2)
  return r2r1
end

push!(fnnames,"FindMach|Find Mach Number from Parameter Ratios|")
"""

### Syntax
    FindMach(which,gam,input)
    FindMach(which,input)

### Description
Calculates the Mach Number of a given shock or flow from the shock jump or stagnation condition relation selected by 'which' for an ideal gas with specific heat ratio 'gam'. Default 'gam' corresponds to Air.

### Arguments
'which::String': Flow Relation Input Selection 
- "PtP" -> Stagnation to static pressure ratio, 'Pt/P' _Eqn. 5.85_
- "TtT" -> Stagnation to static temperature ratio, 'Tt/T' _Eqn. 5.77_
- "M" -> Mach Number immediately after shock, 'M2' _Eqn. 9.60_
- "T" -> Static temperature ratio across shock, 'T2/T1' _Eqn. 9.73_
- "P" -> Static pressure ratio across shock, 'P2/P1' _Eqn. 9.56_
- "Pt" -> Stagnation pressure ratio across shock, 'Pt2/Pt1' _Eqn. 9.66_
- "r" -> Density ratio across shock, 'r2/r1' _Eqn. 9.97_
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Shock Mach Number (Mach Number immediately before shock)

### See Also
PtP, TtT, NormShock

"""
function FindMach(which,gamma,input)
  Mguess = (1+sqrt(5))/2
  if which == "PtP"
    return (2/(gamma-1)*(input^((gamma-1)/gamma) - 1))^(1/2) #Eq 5.85
  elseif which == "TtT"
    return (2/(gamma-1)*(input - 1))^(1/2) #Eq 5.77
  elseif which == "P"
    return (1/(2*gamma)*(input*(gamma+1) + gamma - 1))^(1/2) #Eq 9.56
  elseif which == "Pt"
    FindMach(which,gamma,input,Mguess)
  elseif which == "T"
    FindMach(which,gamma,input,Mguess)
  elseif which == "r"
    FindMach(which,gamma,input,Mguess)
  elseif which == "M"
    FindMach(which,gamma,input,Mguess)
  end
end
function FindMach(which,input)
  FindMach(which,1.4,input)
end
function FindMach(which,gamma,input,Mguess)
  if which == "Pt"
    curval = NormShock("Pt",gamma,Mguess)
    if abs(curval - input)/input < 1e-5
      return round(Mguess,digits=4)
    elseif curval < input
      FindMach(which,gamma,input,0.9*Mguess)
    elseif curval > input
      FindMach(which,gamma,input,1.1*Mguess)
    end
  elseif which == "T"
    curval = NormShock("T",gamma,Mguess)
    if abs(curval - input)/input < 1e-5
      return round(Mguess,digits=4)
    elseif curval > input
      FindMach(which,gamma,input,0.9*Mguess)
    elseif curval < input
      FindMach(which,gamma,input,1.1*Mguess)
    end
  elseif which == "r"
    curval = NormShock("r",gamma,Mguess)
    if abs(curval - input)/input < 1e-5
      return round(Mguess,digits=4)
    elseif curval > input
      FindMach(which,gamma,input,0.9*Mguess)
    elseif curval < input
      FindMach(which,gamma,input,1.1*Mguess)
    end
  elseif which == "M"
    curval = NormShock("M",gamma,Mguess)
    if abs(curval - input)/input < 1e-5
      return round(Mguess,digits=4)
    elseif curval < input
      FindMach(which,gamma,input,0.9*Mguess)
    elseif curval > input
      FindMach(which,gamma,input,1.1*Mguess)
    end
  end
end

push!(fnnames,"OblShock|Oblique Shock Jump Relations|(12.7-15)")
"""

### Syntax
    OblShock(which,M,beta)
    OblShock(which,M,theta)
    OblShock(which,gam,M,beta)
    OblShock(which,gam,M,theta)
    OblShock(which,gam,M,beta,theta)

### Description
Calculates the oblique shock jump relation selected by 'which' for an ideal gas with specific heat ratio 'gam' at shock Mach Number 'M', shock angle 'beta', and flow turn angle 'theta'. Default 'gam' corresponds to Air.

### Arguments
'which::String': Shock Jump Relation Selection 
- "theta"/"θ" -> Flow turn angle through shock, 'theta' _Eqn. 12.15_
- "T" -> Static temperature ratio across shock, 'T2/T1' _Eqn. 12.7_
- "P" -> Static pressure ratio across shock, 'P2/P1' _Eqn. 12.7_
- "Pt" -> Stagnation pressure ratio across shock, 'Pt2/Pt1' _Eqn. 12.7_
- "r" -> Density ratio across shock, 'r2/r1' _Eqn. 12.7_
- "betaS"/"βS" -> Strong shock solution for shock angle, 'beta' _Eqn. 12.15_
  - Requires Syntax: OblShock(which,gam,M,theta)
- "betaW"/"βW" -> Weak shock solution for shock angle, 'beta' _Eqn. 12.15_
  - Requires Syntax: OblShock(which,gam,M,theta)
- "M" -> Mach number immediately after shock, 'M2' _Eqn. 12.7_
  - Requires Syntax: OblShock(which,gam,M,beta,theta)
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta::Real': Shock angle to incoming flow, degrees

'theta::Real': Flow turn angle through shock, degrees

### See Also
NormShock, Mach2Obl, T2T1Obl, P2P1Obl, Pt2Pt1Obl, r2r1Obl, deflang, deflangsol

"""
function OblShock(which::String,gamma::Real,Mach::Real,beta::Real)
  theta = beta
  if which == "theta" || which == "θ"
    deflang(gamma,Mach,beta)
  elseif which == "T"
    T2T1Obl(gamma,Mach,beta)
  elseif which == "P"
    P2P1Obl(gamma,Mach,beta)
  elseif which == "Pt"
    Pt2Pt1Obl(gamma,Mach,beta)
  elseif which == "r"
    r2r1Obl(gamma,Mach,beta)
  elseif which == "betaS" || which == "βS"
    deflangsol("strong",gamma,Mach,theta)
  elseif which == "betaW" || which == "βW"
    deflangsol("weak",gamma,Mach,theta)
  end
end
function OblShock(which::String,gamma::Real,Mach::Real,beta::Real,theta::Real)
  if which == "M"
    Mach2Obl(gamma,Mach,beta,theta)
  end
end
function OblShock(which::String,Mach::Real,beta::Real)
  OblShock(which,1.4,Mach,beta)
end

push!(fnnames,"deflang|**Oblique Shock Flow Turn Angle θ (deg)|(12.15)")
"""

### Syntax
    deflang(gam,M,beta)

### Description
_Eqn. 12.15_

Calculates the flow turn angle through an oblique shock of Mach Number 'M' and angle 'beta' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta::Real': Shock angle to incoming flow, degrees

### See Also
OblShock, deflangsol, Mach2Obl, T2T1Obl, P2P1Obl, Pt2Pt1Obl, r2r1Obl

"""
function deflang(gamma::Real,Mach::Real,shockang::Real)
  numer = Mach^2*sind(shockang)^2 - 1
  denom = 1 + (gamma+1)/2*Mach^2 - Mach^2*sind(shockang)^2
  rhs = cotd(shockang)*(numer/denom)
  theta = atand(rhs)
  return theta
end
#for wolfram alpha atand(cotd(b)*((M^2*sind(b)^2 - 1)/(1 + (gam+1)/2*M^2 - M^2*sind(b)^2)))

push!(fnnames,"deflangsol|**Oblique Shock Angle Solver|(12.15)")
"""

### Syntax
    deflangsol(which,gam,M,beta_guess,theta)
    deflangsol(which,gam,M,theta)

### Description
_Eqn. 12.15_

Calculates the angle of an oblique shock with a flow turn angle 'theta' relative to an incoming flow of Mach Number 'M' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'which::String': Weak or strong shock solution
- "strong" -> strong shock solution
- "weak" -> weak shock solution
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta_guess::Real': Guess for shock angle to incoming flow, degrees

'theta::Real': Flow turn angle through shock, degrees

### See Also
OblShock, deflang, Mach2Obl, T2T1Obl, P2P1Obl, Pt2Pt1Obl, r2r1Obl

"""
function deflangsol(which::String,gamma::Real,Mach::Real,bguess::Real,knownang::Real)
  curval = deflang(gamma,Mach,bguess)
  if which == "weak"
    if abs(curval - knownang)/knownang < 1e-4
      return round(bguess,digits=2)
    elseif curval > knownang
      deflangsol(which,gamma,Mach,0.9*bguess,knownang)
    elseif curval < knownang
      deflangsol(which,gamma,Mach,1.1*bguess,knownang)
    end
  elseif which == "strong"
    if abs(curval - knownang)/knownang < 1e-4
      return round(bguess,digits=2)
    elseif curval < knownang
      deflangsol(which,gamma,Mach,0.9*bguess,knownang)
    elseif curval > knownang
      deflangsol(which,gamma,Mach,1.1*bguess,knownang)
    end
  end
end
function deflangsol(which::String,gamma::Real,Mach::Real,knownang::Real) #Alternative Method to dircumvent initial guess
  if which == "weak"
    deflangsol(which,gamma,Mach,10,knownang)
  elseif which == "strong"
    deflangsol(which,gamma,Mach,90,knownang)
  end
end

push!(fnnames,"Mach2Obl|**Oblique Shock Exit Mach|(12.7)")
"""

### Syntax
    Mach2Obl(gam,M,beta,theta)

### Description
_Eqn. 12.7_

Calculates the Mach Number immediately after an oblique shock of Mach Number 'M', shock angle 'beta', and flow turn angle 'theta' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta::Real': Shock angle to incoming flow, degrees

'theta::Real': Flow turn angle through shock, degrees

### See Also
OblShock, Mach2Norm, deflang, deflangsol, T2T1Obl, P2P1Obl, Pt2Pt1Obl, r2r1Obl

"""
function Mach2Obl(gamma::Real,Mach::Real,beta::Real,theta::Real)
  numer = (gamma-1)*Mach^2*sind(beta)^2 + 2
  denom = 2*gamma*Mach^2*sind(beta)^2 - (gamma-1)
  frac = numer/denom
  M2 = sqrt(frac/sind(beta-theta)^2)
  return M2
end

push!(fnnames,"T2T1Obl|**Oblique Shock Static Temperature Ratio|(12.7)")
"""

### Syntax
    T2T1Obl(gam,M,beta)

### Description
_Eqn. 12.7_

Calculates the static temperature ratio across an oblique shock of Mach Number 'M' and shock angle 'beta' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta::Real': Shock angle to incoming flow, degrees

### See Also
OblShock, T2T1Norm, deflang, deflangsol, Mach2Obl, P2P1Obl, Pt2Pt1Obl, r2r1Obl

"""
function T2T1Obl(gamma::Real,Mach::Real,beta::Real)
  msin = Mach^2*sind(beta)^2
  numer = (2*gamma*msin - (gamma-1)) * ((gamma-1)*msin + 2)
  denom = (gamma+1)^2*msin
  T2T1 = numer/denom
  return T2T1
end

push!(fnnames,"P2P1Obl|**Oblique Shock Static Pressure Ratio|(12.7)")
"""

### Syntax
    P2P1Obl(gam,M,beta)

### Description
_Eqn. 12.7_

Calculates the static pressure ratio across an oblique shock of Mach Number 'M' and shock angle 'beta' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta::Real': Shock angle to incoming flow, degrees

### See Also
OblShock, P2P1Norm, deflang, deflangsol, Mach2Obl, T2T1Obl, Pt2Pt1Obl, r2r1Obl

"""
function P2P1Obl(gamma::Real,Mach::Real,beta::Real)
  numer = 2*gamma*Mach^2*sind(beta)^2 - (gamma-1)
  denom = gamma+1
  P2P1 = numer/denom
  return P2P1
end

push!(fnnames,"Pt2Pt1Obl|**Oblique Shock Stagnation Pressure Ratio|(12.8)")
"""

### Syntax
    Pt2Pt1Obl(gam,M,beta)

### Description
_Eqn. 12.8_

Calculates the stagnation pressure ratio across an oblique shock of Mach Number 'M' and shock angle 'beta' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta::Real': Shock angle to incoming flow, degrees

### See Also
OblShock, Pt2Pt1Norm, deflang, deflangsol, Mach2Obl, T2T1Obl, P2P1Obl, r2r1Obl

"""
function Pt2Pt1Obl(gamma::Real,Mach::Real,beta::Real)
  numer1 = (gamma+1)*Mach^2*sind(beta)^2
  denom1 = (gamma-1)*Mach^2*sind(beta)^2 + 2
  part1 = (numer1/denom1)^(gamma/(gamma-1))
  numer2 = gamma+1
  denom2 = 2*gamma*Mach^2*sind(beta)^2 - (gamma-1)
  part2 = (numer2/denom2)^(1/(gamma-1))
  Pt2Pt1 = part1*part2
  return Pt2Pt1
end

push!(fnnames,"r2r1Obl|**Oblique Shock Density Ratio|(12.7)")
"""

### Syntax
    r2r1Obl(gam,M,beta)

### Description
_Eqn. 12.7_

Calculates the density ratio across an oblique shock of Mach Number 'M' and shock angle 'beta' in an ideal gas with specific heat ratio 'gam'.

### Arguments
'gam::Real': Specific Heat Ratio

'M::Real': Shock Mach Number (Mach Number immediately before shock)

'beta::Real': Shock angle to incoming flow, degrees

### See Also
OblShock, r2r1Norm, deflang, deflangsol, Mach2Obl, T2T1Obl, P2P1Obl, Pt2Pt1Obl

"""
function r2r1Obl(gamma::Real,Mach::Real,beta::Real)
  numer = (gamma+1)*Mach^2*sind(beta)^2
  denom = (gamma-1)*Mach^2*sind(beta)^2 + 2
  r2r1 = numer/denom
  return r2r1
end

push!(fnnames,"aratio|A*/A = f(M)|(10.16)")
"""

### Syntax
    aratio(gam,M)
    aratio(M)
    aratio(which,gam,A*/A)
    aratio(which,A*/A)

### Description
_Eqn. 10.16_

Also known as f(M). Calculates the area ratio between a theoretical choke point and a point at Mach Number 'M' in a flow of an ideal gas with specific heat ratio 'gam'. When used with 'which' and 'A*/A' inputs, will return the Mach Number of the given flow. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Mach Number

'which::String':
- "sup" -> Supersonic Mach Number solution for given A*/A
- "sub" -> Subsoic Mach Number solution for given A*/A

'A*/A::Real': Area ratio or f(M) from which to derive a Mach Number

### See Also
aratiosol, Fanno, Rayleigh, spdofsnd, NormShock, OblShock

"""
function aratio(gamma::Real,Mach::Real) #eq 10.16
  expo = (gamma+1)/(2*(gamma-1))
  numer = Mach*((gamma+1)/2)^expo
  denom = (1 + (gamma-1)/2*Mach^2)^expo
  aratio = numer/denom
  return aratio
end
function aratio(Mach::Real) #Alternative Method for generic case (Air/diatomic)
  aratio(1.4,Mach)
end
# push!(fnnames,"aratio(which,A*/A)|A*/A = f(M) Solver|(10.16)")
function aratio(which::String,gamma::Real,AstarA::Real) #solver
  aratiosol(which,gamma,AstarA,(1+sqrt(5))/2) #arbitrary guess value
end
function aratio(which::String,AstarA::Real) #Alternative Solver Method for generic case (Air/diatomic)
  aratiosol(which,1.4,AstarA,(1+sqrt(5))/2)
end

push!(fnnames,"aratiosol|A*/A = f(M) Solver (included in aratio)|(10.16)")
"""

### Syntax
    aratiosol(which,gam,A*/A,Mguess)

### Description
_Eqn. 10.16_

Calculates the Mach Number of a flow of an ideal gas with specific heat ratio 'gam' from the A*/A, or f(M) of the given flow.

### Arguments
'which::String':
- "sup" -> Supersonic Mach Number solution for given A*/A
- "sub" -> Subsoic Mach Number solution for given A*/A
'gam::Real': Specific Heat Ratio

'A*/A::Real': Area ratio or f(M) from which to derive a Mach Number

'Mguess::Real': Guess for Mach Number (does not affect solution)

### See Also
aratio, Fanno, Rayleigh, spdofsnd, NormShock, OblShock

"""
function aratiosol(which::String,gamma::Real,AstarA::Real,Mguess::Real)
  global curval = aratio(gamma,Mguess)
  if which == "sup"
    if abs(curval - AstarA)/AstarA < 1e-4
      return round(Mguess,digits=4)
    elseif curval < AstarA
      aratiosol(which,gamma,AstarA,0.9*Mguess)
    elseif curval > AstarA
      aratiosol(which,gamma,AstarA,1.1*Mguess)
    end
  elseif which == "sub"
    if abs(curval - AstarA)/AstarA < 1e-4
      return round(Mguess,digits=4)
    elseif curval > AstarA
      aratiosol(which,gamma,AstarA,0.9*Mguess)
    elseif curval < AstarA
      aratiosol(which,gamma,AstarA,1.1*Mguess)
    end
  end
end

push!(fnnames,"Fanno|Fanno Line Relations|(11.31)")
"""

### Syntax
    Fanno(which,gam,M)
    Fanno(which,M)

### Description
Calculates the Fanno Line Flow relation (1D, adiabatic, constant area flow with friction) selected by 'which' for an ideal gas with specific heat ratio 'gam' at Mach Number 'M'. Default 'gam' corresponds to Air.

### Arguments
'which::String': Fanno Line Relation Selection 
- "U" -> Velocity ratio, 'U/U*' _Eqn. 11.31_
- "T" -> Static temperature ratio, 'T*/T' _Eqn. 11.31_
- "P" -> Static pressure ratio, 'P*/P' _Eqn. 11.31_
- "Pt" -> Stagnation pressure ratio, 'Pt*/Pt' _Eqn. 11.31_
- "r" -> Density ratio, 'r*/r' _Eqn. 11.31_
- "CfLD" -> '4CfLmax/D' _Eqn. 11.28_
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Mach Number

### See Also
Rayleigh, aratio, NormShock, OblShock

"""
function Fanno(ver::String,gamma::Real,Mach::Real) #eq 11.31
  part1 = 1 + (gamma-1)/2*Mach^2
  part2 = (gamma+1)/2
  base = part1/part2

  if ver == "U" || ver == "r"
    out = sqrt(Mach^2/base)
  elseif ver == "T"
    out = base
  elseif ver == "P"
    out = Mach*sqrt(base)
  elseif ver == "Pt"
    out = aratio(gamma,Mach)
  elseif ver == "CfLD"
    out = ((1 - Mach^2)/(gamma*Mach^2)) + ((gamma+1)/(2*gamma)) * log(Mach^2/base)
  end
  return out
end
function Fanno(ver::String,Mach::Real) #Alternative Method for generic case (Air/diatomic)
  Fanno(ver,1.4,Mach)
end

push!(fnnames,"Rayleigh|Rayleigh Line Relations|(11.47)")
"""

### Syntax
    Rayleigh(which,gam,M)
    Rayleigh(which,M)

### Description
Calculates the Rayleigh Line Flow relation (1D, frictionless, constant area flow with heat transfer) selected by 'which' for an ideal gas with specific heat ratio 'gam' at Mach Number 'M'. Default 'gam' corresponds to Air.

### Arguments
'which::String': Rayleigh Line Relation Selection 
- "Tt" -> Stagnation temperature ratio, 'Tt*/Tt' _Eqn. 11.43_
- "U" -> Velocity ratio, 'U/U*' _Eqn. 11.47_
- "T" -> Static temperature ratio, 'T*/T' _Eqn. 11.47_
- "P" -> Static pressure ratio, 'P*/P' _Eqn. 11.47_
- "Pt" -> Stagnation pressure ratio, 'Pt*/Pt' _Eqn. 11.47_
- "r" -> Density ratio, 'r*/r' _Eqn. 11.47_
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Mach Number

### See Also
Fanno, aratio, NormShock, OblShock

"""
function Rayleigh(ver::String,gamma::Real,Mach::Real) #eq 11.47
  part1 = 1 + (gamma-1)/2*Mach^2
  base = (1+gamma*Mach^2)/(gamma+1)

  if ver == "U" || ver == "r"
    out = Mach^2/base
  elseif ver == "T"
    out = (base/Mach)^2
  elseif ver == "P"
    out = base
  elseif ver == "Pt"
    out = base*((gamma+1)/(2*part1))^(gamma/(gamma-1))
  elseif ver == "Tt"
    out = (base*(1+gamma*Mach^2))/(2*Mach^2*part1)
  end
  return out
end
function Rayleigh(ver::String,Mach::Real) #Alternative Method for generic case (Air/diatomic)
  Rayleigh(ver,1.4,Mach)
end

push!(fnnames,"PMang|Prandtl-Meyer Expansion Angle ω (deg)|(12.38)")
"""

### Syntax
    PMang(gam,M)
    PMang(M)

### Description
_Eqn. 12.38_

Calculates the Prandtl-Meyer expansion angle 'ω' (degrees) required to reach Mach Number 'M' in a flow of an ideal gas with specific heat ratio 'gam'. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'M::Real': Mach Number

### See Also
PMangsol, aratio, Fanno, Rayleigh, spdofsnd, NormShock, OblShock

"""
function PMang(gamma::Real,Mach::Real)
  frac = ((gamma+1)/(gamma-1))
  part1 = sqrt(Mach^2 - 1)
  PM = sqrt(frac) * atan(sqrt(1/frac)*part1) - atan(part1)
  return PM*180/pi
end
function PMang(Mach::Real) #Alternative Method for generic case (Air/Diatomic)
  PMang(1.4,Mach)
end
#for wolfram alpha sqrt((1.4+1)/(1.4-1)) * atan(sqrt((1.4-1)/(1.4+1))*sqrt(M^2 - 1)) - atan(sqrt(M^2 - 1))

push!(fnnames,"PMangsol|Prandtl-Meyer Expansion Angle Solver|(12.38)")
"""

### Syntax
    PMangsol(gam,Mguess,ω)
    PMangsol(gam,ω)
    PMangsol(ω)

### Description
_Eqn. 12.38_

Calculates the Mach Number reached by expanding a Mach 1 flow through the Prandtl-Meyer expansion angle 'ω' (degrees) in a flow of an ideal gas with specific heat ratio 'gam'. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'Mguess::Real': Mach Number guess

'ω::Real': Prandtl-Meyer Expansion Angle, degrees

### See Also
PMang, aratio, Fanno, Rayleigh, spdofsnd, NormShock, OblShock

"""
function PMangsol(gamma::Real,Mguess::Real,knownang::Real)
  curval = PMang(gamma,Mguess)
  if abs(curval - knownang)/knownang < 1e-5
    return round(Mguess,digits=4)
  elseif curval > knownang
    PMangsol(gamma,0.9*Mguess,knownang)
  elseif curval < knownang
    PMangsol(gamma,1.1*Mguess,knownang)
  end
end
function PMangsol(gamma::Real,knownang::Real) #Alternative Method to circumvent guess
  PMangsol(gamma,10,knownang)
end
function PMangsol(knownang::Real) #Alternative Method for generic case (Air/Diatomic)
  PMangsol(1.4,10,knownang)
end

push!(fnnames,"ShockTubeP4P1|Shock Tube Initial Pressure Ratio|(13.64)")
"""

### Syntax
    ShockTubeP4P1(gam1,gam4,T1,T4,Mw1,Mw4,P2/P1)

### Description
_Eqn. 13.64_

Calculates the initial pressure ratio between the driver and driven sections of a shock tube from the specific heat ratios, molecular weights, and initial static temperatures of both gases, and the pressure ratio ("strength") of the initial shock.

### Arguments
'gam1::Real': Specific Heat Ratio of driven gas

'gam4::Real': Specific Heat Ratio of driver gas

'T1::Real': Initial static temperature of driven gas, K

'T4::Real': Initial static temperature of driver gas, K

'Mw1::Real': Molecular weight of driven gas, amu or kg/mol

'Mw4::Real': Molecular weight of driver gas, amu or kg/mol

'P2/P1::Real': Static pressure ratio across initial shock (shock strength)

### See Also
NormShock, spdofsnd

"""
function ShockTubeP4P1(gam1::Real,gam4::Real,T1::Real,T4::Real,Mw1::Real,Mw4::Real,P2P1::Real)
  numer = (gam4-1)*(spdofsnd(gam1,Mw1,T1)/spdofsnd(gam4,Mw4,T4))*(P2P1-1)
  denom = (4*gam1^2+2*gam1*(gam1+1)*(P2P1-1))^(1/2)
  expo = -(2*gam4/(gam4-1))
  P4P1 = P2P1*(1-numer/denom)^expo
  return P4P1
end

push!(fnnames,"ShockTubeUp|Shock Tube Piston Velocity Up|(13.59)")
"""

### Syntax
    ShockTubeUp(gam1,M1,Mw1,T1)

### Description
_Eqn. 13.59_

Calculates the piston (contact surface) velocity of a shock tube from the specific heat ratio, molecular weight, and initial static temperature of the driven gas, and the Mach Number of the initial shock.

### Arguments
'gam1::Real': Specific Heat Ratio of driven gas

'M1::Real': Mach Number of initial shock

'Mw1::Real': Molecular weight of driven gas, amu or kg/mol

'T1::Real': Initial static temperature of driven gas, K

### See Also
NormShock, spdofsnd

"""
function ShockTubeUp(gamma::Real,Mach::Real,Mw1::Real,Temp::Real)
  Up = spdofsnd(gamma,Mw1,Temp)*((Mach^2-1)/((gamma+1)/2*Mach))
  return Up
end

# push!(fnnames,"ShockTubeMrs(gam1,Up,Mw1,T2,Mguess)|Shock Tube Reflected Shock Mach Number|(13.59)")
# function ShockTubeMrs(gamma::Real,Up::Real,Mw1::Real,Temp::Real,Mguess::Real)
#   curval = ShockTubeUp(gamma,Mguess,Mw1,Temp)
#   if abs(curval - Up)/Up < 1e-5
#     return round(Mguess,digits=4)
#   elseif curval > Up
#     ShockTubeMrs(gamma,Up,Mw1,Temp,0.9*Mguess)
#   elseif curval < Up
#     ShockTubeMrs(gamma,Up,Mw1,Temp,1.1*Mguess)
#   end
# end

push!(fnnames,"ChokedCond|Pressure Ratio for Choked Flow Pt/Pa|(10.26)")
"""

### Syntax
    ChokedCond(gam)
    ChokedCond()

### Description
_Eqn. 10.26_

Calculates the ratio of upstream stagnation pressure to ambient static pressure across a throat required to create choked flow in an ideal gas with specific heat ratio 'gam'. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

### See Also
NozzleMach, gamma, PtP

"""
function ChokedCond(gamma::Real)
  PtPa = ((gamma+1)/2)^(gamma/(gamma-1))
  return PtPa
end
function ChokedCond() #Alternative Method for generic case (Air/diatomic)
  ChokedCond(1.4)
end

push!(fnnames,"NozzleMach|Fully Expanded Nozzle Exit Mach Number|(10.25)")
"""

### Syntax
    ChokedCond(gam,Pt/Pa)
    ChokedCond(Pt/Pa)

### Description
_Eqn. 10.25_

Calculates the Mach Number at the exit of a fully expanded (Pe = Po) converging-diverging nozzle discharging an ideal gas with specific heat ratio 'gam' from the ratio of upstream stagnation pressure to ambient static pressure 'Pt/Pa'. Default 'gam' corresponds to Air.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'Pt/Pa::Real': Ratio of upstream stagnation pressure to ambient static pressure

### See Also
ChokedCond, gamma, PtP

"""
function NozzleMach(gamma::Real,PtPa::Real)
  Me = ((2/(gamma-1)) * (PtPa^((gamma-1)/gamma) - 1))^(1/2)
  return Me
end
function NozzleMach(PtPa::Real) #Alternative Method for generic case (Air/diatomic)
  NozzleMach(1.4,PtPa)
end

push!(fnnames,"scaleht|Atmospheric Scale Height, m|(2.119)")
"""

### Syntax
    scaleht(Mw,T)
    scaleht()

### Description
_Eqn. 2.119_

Calculates the scale height of an atmosphere composed of an ideal gas with molecular weight 'Mw' and atmosphere base temperature 'T'. Defaults correspond to the Earth's atmosphere.

### Arguments
'gam::Real': Specific Heat Ratio (Default: 7/5)

'T::Real': Static temperature at base of atmosphere (Default: 292.15)

### See Also
spdofsnd, gamma

"""
function scaleht(Mw::Real,Temp::Real) #Scale Height of atmoshpere (eq 2.119)
  Ru = 8314.46 #Universal Gas Constant, J/kmol K
  g0 = 9.80665 #Std. Gravitational Acceleration
  H = (Ru*Temp)/(Mw*g0)
  return H
end
function scaleht()
  scaleht(28.97,292.15) #Alternative Method for generic case (Earth)
end

push!(fnnames,"shockRe|Shock Reynolds Number (ρ*U*δ)/μ|(9.96)")
"""

### Syntax
    shockRe(gam,M1,Pr,μᵥ/μ)

### Description
_Eqn. 9.96_

Calculates the shock Reynolds number, (ρ*U*δ)/μ, of a shock of Mach Number 'M' where its thickness δ is the characteristic length, in terms of the specific heat ratio 'gam', Prandtl Number 'Pr', and ratio of bulk to shear (dynamic) viscosity 'μᵥ/μ' of the gas in which it propagates. Shock thickness can be estimated from the shock Reynolds number using the average values of ρ, U, and μ through the shock.

### Arguments
'gam::Real': Specific Heat Ratio

'M1::Real': Shock Mach Number (Mach Number immediately before shock)

'Pr::Real': Prandtl Number

'μᵥ/μ::Real': Ratio of bulk to shear (dynamic) viscosity

### See Also
NormShock, spdofsnd

"""
function shockRe(gamma::Real,Mach::Real,Pr::Real,muvr::Real) #eq 9.96
  numer1 = 2*gamma*(gamma-1)*(4.0/3.0+muvr)*Mach^2*(1/(r2r1Norm(gamma,Mach))-1)^2
  numer2 = (4*(gamma/Pr)*(T2T1Norm(gamma,Mach)-1)^2)/(T2T1Norm(gamma,Mach)+1)
  denom1 = (T2T1Norm(gamma,Mach) + 1)
  denom2 = log(T2T1Norm(gamma,Mach)*(1/r2r1Norm(gamma,Mach))^(gamma-1))
  Re = (numer1 + numer2)/(denom1*denom2)
  return Re
end

push!(fnnames,"CJDetMach|Chapman-Jouget Detonation Mach Number|(11.61)")
"""

### Syntax
    CJDetMach(gam,Δht,Cp,T1)

### Description
_Eqn. 11.61_

Calculates the Mach Number of the shock at the front of a Chapman-Jouget detonation wave from the stagnation enthalpy change produced by the detonation wave 'Δht' and the specific heat ratio 'gam', the specific heat 'Cp', and the static temperature 'T1' of the ambient gas.

### Arguments
'gam::Real': Specific Heat Ratio of ambient gas

'Δht::Real': Stagnation enthalpy change through detonation wave

'Cp::Real': Specific heat capacity of ambient gas

'T1::Real': Static temperature of ambient gas

### See Also
NormShock, spdofsnd, gamma, C_p

"""
function CJDetMach(gamma::Real,dht::Real,Cp_::Real,Temp::Real)
  part1 = 1+(gamma+1)*dht/(Cp_*Temp)
  M1 = (part1 + (part1^2 - 1)^(1/2))^(1/2)
  return M1
end


#~~~~~~~~~~Assemble Documentation From Script Content~~~~~~~~~~
fnnames = [split(fnnames[i],"|") for i=1:length(fnnames)]

local maxl = [maximum([length(fnnames[j][i]) for j=1:length(fnnames)]) for i=1:3]

insert!(fnnames,1,[lpad("",maxl[i],"-") for i=1:3])
insert!(fnnames,3,[lpad("",maxl[i],"-") for i=1:3])
push!(fnnames,[lpad("",maxl[i],"-") for i=1:3])

global aa210docs = [] #for on-demand function reference

for i=1:length(fnnames)
  for j=1:3
    if j==1
      fnnames[i][j] = lpad(fnnames[i][j],maxl[j])
    else
      fnnames[i][j] = rpad(fnnames[i][j],maxl[j])
    end
  end
  push!(aa210docs," "*join(fnnames[i]," | "))
end

# push!(aa210docs, "\n Function Arguments:")

# local pad = 9 #padding width for arg names

#add arg help lines to docstring
# push!(aa210docs, "\n"*lpad("dof: ",pad)*"degrees of freedom of gas molecule (monatomic: 3, diatomic: 5)")
# push!(aa210docs, lpad("Mw: ",pad)*"molar mass of gas, g/mol")
# push!(aa210docs, lpad("gam: ",pad)*"γ, specific heat ratio of gas")
# push!(aa210docs, lpad("T: ",pad)*"temperature of gas (atmosphere base temp for scale height)")
# push!(aa210docs, lpad("M (M1): ",pad)*"mach number of flow (before shock)")
# push!(aa210docs, lpad("Pr: ",pad)*"Prandtl Number of gas")
# push!(aa210docs, lpad("μᵥ/μ: ",pad)*"ratio of bulk viscosity to dynamic viscosity")

# push!(aa210docs, lpad("which: ",pad)*"flow relation selection")
# push!(aa210docs, lpad(" ",pad+2)*"\"U\" -> U/U*, ρ*/ρ (Rayleigh & Fanno Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"M\" -> Shock Exit Mach Number M2 (NormShock & OblShock Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"T\" -> T*/T (Rayleigh & Fanno), T2/T1 (NormShock & OblShock)")
# push!(aa210docs, lpad(" ",pad+2)*"\"P\" -> P*/P (Rayleigh & Fanno), P2/P1 (NormShock & OblShock)")
# push!(aa210docs, lpad(" ",pad+2)*"\"Pt\" -> Pt*/Pt (Rayleigh & Fanno), Pt2/Pt1 (NormShock & OblShock)")
# push!(aa210docs, lpad(" ",pad+2)*"\"CfLD\" -> 4CfLmax/D (Fanno only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"Tt\" -> Tt*/Tt (Rayleigh Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"r\" -> ρ2/ρ1 (NormShock & OblShock Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"sub\" -> subsonic solution (aratiosol Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"sup\" -> supersonic solution (aratiosol Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"strong\" -> strong shock solution (deflangsol Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"weak\" -> weak/normal shock solution (deflangsol Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"thetaS\"/\"θS\" -> strong shock solution (OblShock Only)")
# push!(aa210docs, lpad(" ",pad+2)*"\"thetaW\"/\"θW\" -> weak/normal shock solution (OblShock Only)")

# push!(aa210docs, lpad("beta: ",pad)*"oblique shock angle to oncoming flow (degrees)")
# push!(aa210docs, lpad("theta: ",pad)*"flow turn angle thru oblique shock (degrees)")
# push!(aa210docs, lpad("Δht: ",pad)*"change in stagnation enthalpy through detonation")

push!(aa210docs, "\nFunctions with descriptions preceded by an asterisk are included within NormShock")

push!(aa210docs, "\nFunctions with descriptions preceded by two asterisks are included within OblShock")

push!(aa210docs, "\nType \"AA210AFnList()\" or \"?AA210A_Functions\" to see function list again")

push!(aa210docs, "\nType \"?<function_name>\" to see individual function documentation")

insert!(aa210docs,1,"\n\n AA210A Compressible Flow Formulae\n\n Material from \"AA210A Course Reader - Fundamentals of Compressible Flow\"\n by Prof. Brian J. Cantwell - https://web.stanford.edu/~cantwell/\n")

end #close "begin" local variable scope

aa210docscomp = join(aa210docs,"\n") #compile docstring from array

"""

Call AA210AFnList() for function list

"""
function AA210AFnList() #function to print docstring on demand
  println(aa210docscomp)
end

end #module end

#Created by Jeff Robinson - jbrobin@stanford.edu