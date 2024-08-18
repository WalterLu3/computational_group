set t 'stages'    /t1*t3/
    s 'scenarios' /s1*s3/
    j 'iteration index' / j1*j20 /
    i 'trial index' /i1*i2/;

Set
   tt(t) 'control set for stages'
   jj(j) 'dynamic j'
;


Parameter
    inflow(t)
    cost(t)    /t1 50
                t2 100
                t3 150/
;

Parameter
    sto_inflow(s) /s1 0
                   s2 50
                   s3 100/;

Positive Variable
    volume(t)
    thermalGeneration(t)
    hydroGeneration(t)
    hydroSpill(t)
;

Free     Variable ACOST    'approximation of cost';
Positive Variable ALPHA(t) 'approximation of future cost function (FCF)';

volume.up(t) = 200;

Equation
    stateTransfer(t) 'State update variable'
    demandSatisfy(t) 'Need to meet demand'
;

stateTransfer(tt(t))..
    volume(t) =e= volume(t-1) + 200$(sameas('t1',t)) - hydroGeneration(t) - hydroSpill(t) + inflow(t);

demandSatisfy(tt(t))..
    hydroGeneration(t) + thermalGeneration(t) =e= 150;




parameter
    cont_m(j,i,t) 'dual variables associated with the mass balance constraint'
    delta(j,i,t)  'RHS of the Benders cuts';

* no cuts at the beginning 
cont_m(j,i,t) = 0;
delta(j,i,t)  = 0;
jj(j)         = no;

Equation
    Obj_Approx
    Cuts(j,i,t);

Obj_Approx..
    ACOST =e= sum(tt(t), cost(t) * thermalGeneration(t))
              + sum(tt(t), ALPHA(t+1));

* no cuts for the leaf node
Cuts(jj,i,t)$(ord(t) < card(t))..
   ALPHA(t+1) - cont_m(jj,i,t+1)*volume(t) =g= delta(jj,i,t);

Model waterSDDP / all /;

Parameter i_dec(i,t) 'trial decisions';

Set jtis(j,t,i,s) 'scenario trial sample';
execseed = 1 + gmillisec(jnow);
* Select some for the first backward recursion;
loop(i, i_dec(i,t) =60 + (ord(i) - 1)/(card(i) - 1)*(200));
loop(s$sameas('s1',s), jtis(j,t,i,s + UniFormInt(0,card(s) - 1))$(not sameas('t1',t)) = yes);


Set revt(t,t) 'backward loop for backward recursion excluding week 1';
loop(t$(ord(t) < card(t)), revt(t,t + (card(t) - 2*ord(t) + 1)) = yes);


* parameter to store solution in each trial
Parameter
    is_volume(i,s,t)
    is_inflow(i,s,t)
    is_stateTransfer(i,s,t)
    is_demandSatisfy(i,s,t)
    is_Cuts(i,s,j,i,t)
    is_ACOST(i,s)
    is_ALPHA(i,s,t)
    so / SkipBaseCase 1, LogOption 1, UpdateType 2 /;

inflow(t) = 50;

Set
   dict_b 'backward recursion'
          / is             .scenario.''
            so             .opt     .''
            volume         .fixed   .is_volume
            inflow         .param   .is_inflow
            stateTransfer  .marginal.is_stateTransfer
            demandSatisfy  .marginal.is_demandSatisfy
            Cuts           .marginal.is_Cuts/
   dict_f 'forward simulation'
          / is      .scenario.''
            so      .opt     .''
            volume  .fixed   .is_volume
            inflow  .param   .is_inflow
            ACOST   .level   .is_ACOST
            ALPHA   .level   .is_ALPHA
            volume  .level   .is_volume/;


$if not set timelim $set timelim 2 

Scalar start;
start = jnow;

Parameter
   conv(j,*) 'convergence parameters'
   zt(j,i)   'objective for trial solution';

Set is(i,s);

scalar test;

Alias (t,t1,t2);
Alias (i,ip);

loop(j, 
   is(i,s) = yes;      
   loop((t1,t2)$revt(t1,t2),
      tt(t) = no; tt(t2) = yes;;

      option clear = is_volume, clear = is_inflow;

* what is it doing here? fix the last decision
      is_volume(i,s,t2-1) = i_dec(i,t2-1);
      is_inflow(i,s,t2) = sto_inflow(s);

      solve waterSDDP min ACOST using lp scenario dict_b;

      cont_m(j,i,t2)   =  sum(s, is_stateTransfer(i,s,t2))/card(s);
      delta(j,i,t2-1)  = [sum{s, is_stateTransfer(i,s,t2)*is_inflow(i,s,t2)
                                + is_demandSatisfy(i,s,t2)*150}
                                + sum((s,jj,ip),
                                  is_cuts(i,s,jj,ip,t2)*delta(jj,ip,t2))$(ord(t2) < card(t))]/card(s);
* we want cuts from the stage t+1
      jj(j) = yes;
   );

* forward pass
   loop(t1,
      tt(t) = no; tt(t1) = yes;
      if(sameas('t1',t1),
* first stage is done specially
* fix 0 stage volume
        volume.fx(t1-1) = 200;
        solve waterSDDP min ACOST using lp;
* store the first stage cost as the lower bound
        conv(j,'lo') = ACOST.l;
        zt(j,i)      = ACOST.l - ALPHA.l(t1+1);

* update the policy
        i_dec(i,'t1') = volume.l('t1');
        
      else
* other stage done normally
        is(i,s) = jtis(j,t1,i,s); 
        option clear = is_volume, clear = is_inflow;

* get the state variable value from last stage
        is_volume(is(i,s),t1-1) = i_dec(i,t1-1);
* get the random variable value in a specific scenario
        is_inflow(is(i,s),t1) = sto_inflow(s);
        solve waterSDDP min ACOST using lp scenario dict_f;

        i_dec(i,t1) = sum(is(i,s), is_volume(i,s,t1));
        zt(j,i) = zt(j,i) + sum(is(i,s), is_ACOST(i,s) - is_ALPHA(i,s,t1+1));
      );
   );
);

display conv;
display zt;
