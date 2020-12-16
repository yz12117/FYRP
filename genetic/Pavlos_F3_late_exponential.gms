***********************
* OPTIMIZATION FILE 5 *
***********************

***** Specifying the directories where the root files are present
$set myroot InputForGAMS/Pavlos/F3/late_exponential/
$set condition late_exponential
$set feed F3
$set doc Pavlos
*$set doc SAR_FeedC
*$set doc SAR_FeedAll
*$set doc SAR_FeedAll40
*$set doc WID

***** Specifying the minimum threshold of biomass
$set biom 0.016202 * 0.1

Options
* minimize the size of the LST file
        limrow = 0
        limcol = 0
        solprint = off
        sysout = off
* solution and solver options
*        optCR = 1E-9
*        optCA = 1E-9
        iterlim = 1000000
        decimals = 8
        reslim = 100000
        work = 50000000;

******************************** SET DEFINITIONS ***************************************
*
*       i                       = set of all metabolites
*       j                       = set of all reactions
*       offaeroglucose(j)       = set of reactions that should be turned off under aerobic glucose conditions
*       exchange(j)             = set of reactions defined in the minimal M9 media for uptake
*       essential(j)            = set of essential in silico essential reactions under aerobic condition
*       biom(j)                 = set containing reaction for biomass flux
*       nogene(j)               = set of reactions that do not have any gpr associations
*       excg(j)                 = set of transport reactions
*       blocked(j)              = set of blocked reactions
*       index                   = set to store alternate interventions
*
*****************************************************************************************
Sets

i
$include "%myroot%ModelMETS.txt"

j
$include "%myroot%ModelRXNames.txt"

exchange(j)
$include "%myroot%ModelUptakes.txt"

Biom(j) / 'R5261' /

constraint_flux(j)
/
* Biomass
'R5261'

* IgG
'R1082'

* Glutamine synthase
*'R614'
/

*UpregulateThese(j)
*$include "%myroot%ModelEssentialRxns.txt"

NoGenes(j)
$include "%myroot%ModelNoGenes.txt"

index /1*100/
;

****************************** PARAMETRIC CONSTANTS USED ********************************
*
*       s(i,j)          = Stoichiometry of metabolite i in reaction j
*       rxntype(j)      = Specified whether a reaction is irreversible (0), reversible
*                         (1 - forward, 2 - backward reaction), or exchange reactions (4)
*       LB(j)/UB(j)     = Stores the lower and upper bounds for each reaction j
*       low(j)/high(j)  = Stores the minimum/maximum flux od each reaction
*       epsilon         = A small value
*       x               = counter
*       alt(index,j)    = Stores alternate solutions of reaction removals. For each
*                         index, a value of 1 indicates reaction is removed, and a
*                         value of 0 indicates reaction is not removed
*
*****************************************************************************************

Parameters

s(i,j)
$include "%myroot%ModelS.txt"

rxntype(j)
$include "%myroot%ModelRXNStype.txt"

Reverse(j)
$include "%myroot%ModelBackwardRxns.txt"

*** flux bounds
LB(j) Lower bounds containing the medium and the feed
$include "%myroot%ModelLowerBound_%doc%.txt"

UB(j) Upper bounds containing the medium and the feed
$include "%myroot%ModelUpperBound_%doc%.txt"

maxb, maxi, cb, ci
e_up, e_down

*** matrix to store Interventions
done_up(index,j)
done_ko(index,j)
done_down(index,j)

w1,w2,w3
jlnj
null_par(j), up_par(j), down_par(j)
;

***************************** VARIABLE DEFINITIONS **************************************
*       REAL VARIABLES:
*       v(j)            = Flux of a reaction j (+/-)
*       z               = Objective value (+/-)
*       lambda(i)       = Dual variable for Primal constraints on network stoichiometry
*       wL(j)/wU(j)     = Variable for lineraing product of a dual and binary variales
*       POSITIVE VARIABLES:
*       muL(j)/muU(j)   = Dual variables for Primal constraints on limits of reaction fluxes
*       BIANRY VARIABLES
*       y_null(j)            = Indicates 1 for reaction removal, 0 otherwise
*
******************************************************************************************
Variables
z, b
v(j)
lambda(i)
mu(j)
wL(j), wU(j)
wphi(j), wtheta(j)
wFOR(j), wBACK(j)
;

Negative Variables
biomass
;

Positive variables
phi(j), theta(j)
muL(j), muU(j)
mFor(j),mBack(j)
;

binary variables
y_null(j)
y_low(j)
y_up(j)
y_rev(j)
;

**************** INITIALIZING PARAMETRIC VARIABLES AND SETS ****************************

*** Scalars

scalar counter /1/;
scalar M /100/;
scalar k ;
scalar n /0/;
scalar cutoff /0/;
scalar cutoff1 /0/;


*****************************************************************************************
********************************* Equations definitions *********************************
*****************************************************************************************
*
*       OUTER LEVEL CONSTRAINTS:
*       Obj                             = Maximize flux of target chemical
*       Primal_dual                     = Strong duality cosntraint to equate primal objective equal to dual objective
*       Int_cut                         = Integer cut constraint to find alternate solutions
*       Tot                             = Total number of interventions
*       Lin_muLa/b/c/d(j)               = Constraints to linearize muL(j)*y_null(j) with wL(j)
*       Lin_muUa/b/c/d(j)               = Constraints to linearize muU(j)*y_null(j) with wL(j)
*       PRIMAL CONSTRAINTS:
*       Stoic(i)                        = Stoichiometric Constraint
*       Flux_limL(j)/Flux_limU(j)       = Lower/Upper limits on flux v(j)
*       DUAL CONSTRAINTS:
*       dual1a                          = Dual constraint for primal biomass reaction variable
*       dual1b                          = Dual constraint for primal variables other than bioamss
*

Equations

Obj
Stoic
bio, Flux_Upper, Flux_Lower, Flux_Upper_t, Flux_Lower_t
Con_bio
Flux_limL,Flux_limU
Flux_Upregulate
Flux_Downregulate

Primal_dual
dual1a,dual1b,dual1c

Int_cut
Tot
y_lim

back_lin1
back_lin2
back_lin3
back_lin4
RevRXN
ForBackRxn

UP_outer4, UP_outer5, UP_outer6, UP_outer7
DOWN_outer8, DOWN_outer9, DOWN_outer10, DOWN_outer11
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd

lin_con_1
lin_con_2
;

*****************************************************************************************
******************************* Defining the constraints. *******************************
*****************************************************************************************

bio..                                         b =e= v('R5261');
Flux_Upper(j)..                                                           v(j) =l= UB(j);
Flux_Lower(j)..                                                           -v(j) =l= -LB(j);

Flux_Upper_t(j)..                                                         v(j) =l= (1-null_par(j)) * ( UB(j) - (((UB(j) - LB(j))*e_down)*down_par(j)));
Flux_Lower_t(j)..                                                         -v(j) =l= -(1-null_par(j)) * (LB(j) + (((UB(j) - LB(j))*e_up)*up_par(j)));

Con_bio..                                     -v('R5261') =l= -maxb*cb ;
Obj..                                                                             z =e= v('R1082');
Stoic(i)..                                                                        sum(j, (S(i,j)*v(j))) =e= 0;

Flux_Downregulate(j)..                                            v(j) =l= (UB(j) - ((UB(j) - LB(j))*e_down))*y_low(j) + UB(j)*(1-y_low(j)-y_null(j));
Flux_limL(j)..                                                            -v(j) =l= -LB(j)*(1-y_null(j) );
Flux_limU(j)..                                                            v(j) =l= UB(j)*(1-y_null(j) );
Flux_Upregulate(j)..                                              -v(j) =l= -((LB(j) + ((UB(j) - LB(j))*e_up))*y_up(j) + LB(j)*(1-y_up(j)-y_null(j)));

Tot..                                                                         sum(j, y_null(j) + y_low(j)+ y_up(j)) =g= k;
y_lim(j)..                                                                        y_null(j) + y_low(j) + y_up(j) =l= 1;
Int_cut(index)..                                                      sum(j, done_up(index,j) * y_up(j) + done_ko(index,j) * y_null(j)  + done_down(index,j) * y_low(j)) =l= k-1;

RevRXN(j)..                                                             v(j) =l= (1-y_rev(j)) * UB(j);
ForBackRxn(j)$(Reverse(j) > 1)..                y_rev(j) + y_rev(j)$(ord(j)=Reverse(j)) =l= 1;
back_lin1(j)..                                                  wBACK(j) =l= M*(y_rev(j));
back_lin2(j)..                                                  wBACK(j) =g= -M*(y_rev(j));
back_lin3(j)..                                                  wBACK(j) =l= mBack(j) + M*(1-y_rev(j));
back_lin4(j)..                                                  wBACK(j) =g= mBack(j) - M*(1-y_rev(j));

Primal_dual..           v('R5261') =e=  - sum(j, wtheta(j)*(LB(j) + ((UB(j) - LB(j))*e_up)) + theta(j)*LB(j) - wtheta(j)*LB(j) )
                                                                        + sum(j, wphi(j)*(UB(j) - ((UB(j) - LB(j))*e_down)) + phi(j)*UB(j) - wphi(j)*UB(j) )
                                                                        + sum(j, muU(j)*UB(j) - wU(j)*UB(j))
                                                                        - sum(j, muL(j)*LB(j) - wL(j)*LB(j))
                                                                        + sum(j, mBack(j)*UB(j) - wBACK(j)*UB(j));

Lin_muLa(j)..           wL(j) =l= M*y_null(j);
Lin_muLb(j)..           wL(j) =g= -M*y_null(j);
Lin_muLc(j)..       wL(j) =l= muL(j) + M*(1-y_null(j) );
Lin_muLd(j)..           wL(j) =g= muL(j) - M*(1-y_null(j) );

Lin_muUa(j)..       wU(j) =l= M*y_null(j);
Lin_muUb(j)..       wU(j) =g= -M*y_null(j);
Lin_muUc(j)..       wU(j) =l= muU(j) + M*(1-y_null(j) );
Lin_muUd(j)..       wU(j) =g= muU(j) - M*(1-y_null(j) );

DOWN_outer8(j)..    wphi(j) =l= M*y_low(j);
DOWN_outer9(j)..    wphi(j) =g= -M*y_low(j);
DOWN_outer10(j)..   wphi(j) =l= phi(j) + M*(1-y_low(j) );
DOWN_outer11(j)..   wphi(j) =g= phi(j) - M*(1-y_low(j) );

UP_outer4(j)..      wtheta(j) =l= M*y_up(j);
UP_outer5(j)..      wtheta(j) =g= -M*y_up(j);
UP_outer6(j)..      wtheta(j) =l= theta(j) + M*(1-y_up(j) );
UP_outer7(j)..      wtheta(j) =g= theta(j) - M*(1-y_up(j) );



dual1b(j)$(not biom(j))..       sum(i,lambda(i)*S(i,j)) - theta(j) + phi(j) + muU(j) - muL(j) + mBack(j) =e= 0;
dual1c(j)$(biom(j))..           sum(i, lambda(i)*S(i,j)) - theta(j) + phi(j) + muU(j) - muL(j) + mBack(j) =e= 1;

*Primal_dual..          v('R5190') =e=  - sum(j, wtheta(j)*(LB(j) + ((UB(j) - LB(j))*e_up)) + theta(j)*LB(j) - wtheta(j)*LB(j) )
*                                                                       + sum(j, wphi(j)*(UB(j) - ((UB(j) - LB(j))*e_down)) + phi(j)*UB(j) - wphi(j)*UB(j) )
*                                                                       + sum(j, muU(j)*UB(j) - wU(j)*UB(j))
*                                                                       - sum(j, muL(j)*LB(j) - wL(j)*LB(j))
*                                                                       + sum(j, mFor(j)*UB(j) - wFOR(j)*UB(j))
*                                                                       + sum(j, mBack(j)*UB(j) - wBACK(j)*UB(j));

*dual1b(j)$(not biom(j))..      sum(i,lambda(i)*S(i,j)) - theta(j) + phi(j) + muU(j) - muL(j) + mBack(j) + mFor(j) =e= 0;
*dual1c(j)$(biom(j))..          sum(i, lambda(i)*S(i,j)) - theta(j) + phi(j) + muU(j) - muL(j) + mBack(j) + mFor(j) =e= 1;

*Primal_dual..          v('R5190') =e= sum(j, muU(j)*UB(j) - wU(j)*UB(j)) - sum(j, muL(j)*LB(j) - wL(j)*LB(j));
*Primal_dual..          v('R5190') =e= biomass*maxb*cb + sum(j, wtheta(j)*(LB(j) + ((UB(j) - LB(j))*e_up)) + theta(j)*LB(j) - wtheta(j)*LB(j) )
*                                                                       - sum(j, wphi(j)*(UB(j) - (UB(j) - LB(j))*e_down) + phi(j)*UB(j) - wphi(j)*UB(j) )
*                                                                       + sum(j, muU(j)*UB(j) - wU(j)*UB(j))
*                                                                       - sum(j, muL(j)*LB(j) - wL(j)*LB(j));

*dual1a(j)..            sum(i,lambda(i)*S(i,'E39')) + theta('E39') - phi('E39') + muU('E39') - muL('E39')  =e= 1;
*dual1b(j)$(not biom(j))..      sum(i,lambda(i)*S(i,j)) + theta(j) - phi(j) + muU(j) - muL(j)  =e= 0;
*dual1c(j)..                             biomass*maxb*cb + sum(i, lambda(i)*S(i, 'R5190')) + theta('R5190') - phi('R5190') + muU('R5190') - muL('R5190') =e= 1;
*dual1b(j)$(not biom(j))..      sum(i,lambda(i)*S(i,j)) + muU(j) - muL(j)  =e= 0;
*dual1c(j)..                             sum(i, lambda(i)*S(i, 'R5190')) + muU('R13') - muL('R13') =e= 1;

*****************************************************************************************
************** Defining Upper/Lower bounds and fixing the binary variables **************
*****************************************************************************************

scalar max /100/;

LB(j)$(constraint_flux(j)) = 0;
UB(j)$(constraint_flux(j)) = 10;

v.lo(j)$(not constraint_flux(j) ) = 0;
v.up(j)$(not constraint_flux(j) ) = max;
v.lo(j)$(constraint_flux(j) ) = 0;
v.up(j)$(constraint_flux(j) ) = max;

* Fix the forwards and backward binaries
y_rev.fx(j)$(rxntype(j)=0 or rxntype(j)=3 or rxntype(j)=4) = 0;

* No change in the exchange reactions.
y_up.fx(j)$(exchange(j)) = 0;
y_null.fx(j)$(exchange(j)) = 0;
y_low.fx(j)$(exchange(j)) = 0;

* No change in the important reactions.
y_up.fx(j)$(constraint_flux(j)) = 0;
y_null.fx(j)$(constraint_flux(j)) = 0;
y_low.fx(j)$(constraint_flux(j)) = 0;

* No change in the exchange reactions.
y_up.fx(j)$(rxntype(j)=4) = 0;
y_null.fx(j)$(rxntype(j)=4) = 0;
y_low.fx(j)$(rxntype(j)=4) = 0;

* No change in the important reactions.
*y_up.fx(j)$(not UpregulateThese(j) ) = 0;
*y_low.fx(j)$(UpregulateThese(j) ) = 0;
*y_null.fx(j)$(UpregulateThese(j) ) = 0;

* No change in the blocked reactions.
y_up.fx(j)$(LB(j)=0 and UB(j)=0) = 0;
y_low.fx(j)$(LB(j)=0 and UB(j)=0) = 0;
y_null.fx(j)$(LB(j)=0 and UB(j)=0) = 0;

* No change in the reactions without genes and stuff.
y_up.fx(j)$(NoGenes(j)) = 0;
y_low.fx(j)$(NoGenes(j)) = 0;
y_null.fx(j)$(NoGenes(j)) = 0;

* Disable Up and Down regulations
*y_up.fx(j) = 0;
*y_low.fx(j) = 0;

*****************************************************************************************
********************************** DECLARING THE MODEL **********************************
*****************************************************************************************

Model MaxBio
/
Stoic bio Flux_Upper Flux_Lower
/;

Model MaxIgG
/
Stoic Obj Flux_Upper Flux_Lower
/;

Model TestModel
/
Stoic bio Flux_Upper_t Flux_Lower_t
/;

model bilevel
/

Obj
Stoic
*Con_bio

Flux_limL,Flux_limU
Flux_Upregulate
Flux_Downregulate

Primal_dual
dual1b
dual1c

Int_cut
Tot
y_lim

back_lin1
back_lin2
back_lin3
back_lin4
RevRXN
ForBackRxn

UP_outer4, UP_outer5, UP_outer6, UP_outer7
DOWN_outer8, DOWN_outer9, DOWN_outer10, DOWN_outer11
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd
/
;


*bilevel.optfile = 5;
*bilevel.holdfixed = 1;

*****************************************************************************************
*************** Solving the MILP and priting the solution in a .txt file. ***************
*****************************************************************************************
*** 1. Iterate through total number of interventions starting from a single reaction removal
*** 2. Maximize the product flux and print its value and the interventions
*** 3. Store the intervention combination in alt(index,j)
*** 4. If an alternate solution with non-zero flux of product exists, find it;
***    else, move to higher number of intervention

file forced /%doc%_%feed%_%condition%_large_n.txt/;
put forced;
put '****************************************'//;

solve MaxBio using lp maximizing b;
solve MaxIgG using lp maximizing z;
put 'Max Biomass ='b.l:0:8/;
put 'Max GS ='z.l:0:8/;
maxb = b.l;
maxi = z.l;

*** setting LB wild biomass flux to 10% of theoretical maximum
LB(j)$(biom(j)) = %biom%;

cb = 0.1;
ci = 0;

e_up = 0.7;
e_down = 0.3;

*for (e_up = 0.4 to 1.1 by 0.1,
*for (e_down = 0.2 to 1 by 0.25,

cutoff = 0;
cutoff1 = 0;
jlnj = 0;

*for (k = 10 to 100 by 10,
k = 10;
put "Upregulation coefficient is: "e_up, "Downregulation coefficient is: ",e_down /;

        counter = 1;
        n = 0;

    done_up(index,j) = 0;
    done_ko(index,j) = 0;
    done_down(index,j) = 0;

        put 'Number of Interventions (k) =' k/;
        put '**********************************************************'//;
        while( counter = 1,

                solve bilevel using mip maximizing z;
                jlnj = jlnj + 1;
                n = n + 1;
                if(n eq 1, cutoff = z.l );

                if(bilevel.modelstat eq 1 or bilevel.modelstat eq 8 ,
                        if( z.l <= cutoff1  * 1.01,
                                counter = 0;
                                put "Objective value/better solution not found Cutoff1 = ", cutoff1:0:8, "z = ", z.l:0:8 /;
                                cutoff1 = cutoff;
                        );
                                put 'Total Iterations are' jlnj //;

                                put 'Iteration No. =' n//;

                                        put 'Biomass ='v.l('R5261'):0:8/;
                                        put 'v(IgG) = 'z.l:0:8/;
*                                       put 'v(GS) = 'z.l:0:8/;

                                        put '****************************************'//;
                                        put 'Testing the regulations!' /;
                                        up_par(j) = y_up.l(j);
                                        null_par(j) = y_null.l(j);
                                        down_par(j) = y_low.l(j);
                                        solve TestModel using lp maximizing b;
                                        put 'Biomass is ='b.l:0:8/;
                                        put 'IgG flux is = 'v.l('R1082'):0:8/;
*                                       put 'Glutamine Synthase flux is = 'v.l('R614'):0:8/;
                                        put '****************************************'//;

                                        put /'Gene Regulations:'//;
                                        put "%%%%% Increase the:"/;
                                        w1 = 0;
                                        loop (j$(y_up.l(j) = 1),
                                                w1 = w1 + 1;
                                                put "r",w1:0:0,"=",ord(j):0:0,";"/;
                                        );
                                        put "b = [];"/"for i=1:"w1:0:0/"eval(sprintf('b = [b; r%d];',i));"/"end"/"r_up = b;"/;
                                        w2 = w1;
                                        w1 = w1 + 1;
                                        put "%%%%% Knock-out the:"/;
                                        loop (j$(y_null.l(j) = 1),
                                                w2 = w2 + 1;
                                                put "r",w2:0:0,"=",ord(j):0:0,";"/;
                                        );
                                        put "b = [];"/"for i="w1:0:0":"w2:0:0/"eval(sprintf('b = [b; r%d];',i));"/"end"/"r_null = b;"/;
                                       w3 = w2;
                                       w2 = w2 + 1;
                                       put "%%%%% Decrease the:"/;
                                       loop (j$(y_low.l(j) = 1),
                                               w3 = w3 + 1
                                               put "r",w3:0:0,"=",ord(j):0:0,";"/;
                                       );
                                       put "b = [];"/"for i="w2:0:0":"w3:0:0/"eval(sprintf('b = [b; r%d];',i));"/"end"/"r_down = b;"//;

                        done_up(index,j)$(ord(index) = n and (y_up.l(j) = 1)) = 1;
                        done_ko(index,j)$(ord(index) = n and (y_null.l(j) = 1)) = 1;
                        done_down(index,j)$(ord(index) = n and (y_low.l(j) = 1)) = 1;



                );

                if(n >= 1 or (bilevel.modelstat ne 8 and bilevel.modelstat ne 1),
                                counter = 0;
                put "No feasible solution achieved, modelstat :"bilevel.modelstat, " n = ",n//;
                                cutoff1 = cutoff;
                );

        put '**************************************************'//;

        );

        put /'**********************************************************'/;
*);
* k for loop
*);
* e_down for loop
*);
* e_up for loop


putclose forced;
*****************************************************************************************
*****************
