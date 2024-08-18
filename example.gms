$ontext
    https://www.gams.com/latest/docs/UG_EMP_SP.html#INDEX_EMP_22_stochastic_21_programming
$offtext

Variable  z      "profit";
Positive Variables
          x      "units bought"
          i      "inventory"
          l      "lost sales"
          s      "units sold";

Scalars   c      "purchase costs per unit"         / 30 /
          p      "penalty shortage cost per unit"  /  5 /
          h      "holding cost per leftover unit"  / 10 /
          v      "revenue per unit sold"           / 60 / 
          d      "demand, random parameter"        / 45 /;

Equations profit "profit to be maximized" 
          row1   "demand = UnitsSold + LostSales"
          row2   "inventory = UnitsBought - UnitsSold";

profit..  z =e= v*s - c*x - h*i - p*l;
row1..    d =e= s + l;
row2..    i =e= x - s;

Model nv / all /;

File emp / '%emp.info%' /;
put emp '* problem %gams.i%'/;
$onput
randvar d discrete 0.7 45 0.2 40 0.1 50
stage 2 i l s d
stage 2 Row1 Row2
$offput
putclose emp;

Set scen           "scenarios" / s1*s3 /;
Parameter
     s_d(scen)     "demand realization by scenario"
     s_x(scen)     "units bought by scenario"
     s_s(scen)     "units sold by scenario"
     s_rep(scen,*) "scenario probability"   / #scen.prob 0/;
 
Set dict / scen .scenario.''
           d    .randvar .s_d
           s    .level   .s_s
           x    .level   .s_x
           ''   .opt     .s_rep /;
 
solve nv max z use EMP scenario dict;
display s_d, s_x, s_s, s_rep;