%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Teoría Macrodinámica %%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%    FIEECS - UNI      %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma CRRA. 
// (c) Carlos Rojas Quiroz 

var lab c w r y kap innv z g;
predetermined_variables kap;
varexo e_z e_g;
parameters alpha delta betta sigma kappa eta rho_z rho_g 
z_ss lab_ss r_ss  kap_ss w_ss y_ss c_ss inv_ss g_ss C_Y I_Y G_Y;

alpha  = 0.6500;
delta  = 0.0250;
betta  = 0.9909;
sigma  = 1;
eta    = 1;
rho_z  = 0.7203;
rho_g  = 0.5715;
z_ss   = 1;
G_Y    = 0.15282;
lab_ss = 0.383409634485720;
y_ss   = z_ss^(1/alpha)*((1-alpha)*betta/(1-betta+betta*delta))^((1-alpha)/alpha)*lab_ss;
kappa  = (alpha*1/(1-((1-alpha)*betta*delta/(1-betta+betta*delta))-G_Y)^sigma)*(y_ss^(1-eta*sigma))/(lab_ss^(1+eta));
w_ss   = alpha*y_ss/lab_ss;
kap_ss = (1-alpha)*betta/(1-betta+betta*delta)*y_ss;
inv_ss = delta*kap_ss;
r_ss   = (1-alpha)*y_ss/kap_ss-delta;
c_ss   = ((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)*y_ss;
g_ss   = G_Y*y_ss;
C_Y    = c_ss/y_ss;
I_Y    = inv_ss/y_ss;

model;
exp(w)          = kappa*(exp(lab)^eta)*(exp(c)^sigma);
exp(c)^(-sigma) = betta*exp(c(+1))^(-sigma)*(1+exp(r(+1)));
exp(w)          = alpha*exp(y)/exp(lab);
exp(r)+delta    = (1-alpha)*exp(y)/exp(kap);
exp(y)          = exp(c)+exp(innv)+exp(g);
exp(kap(+1))    = (1-delta)*exp(kap)+exp(innv);
exp(y)          = exp(z)*exp(kap)^(1-alpha)*exp(lab)^alpha;
z               = (1-rho_z)*log(z_ss) + rho_z*z(-1) + e_z;
g               = (1-rho_g)*log(g_ss) + rho_g*g(-1) + e_g;
end;

steady_state_model;
lab =log(lab_ss);
c   =log(c_ss); 
w   =log(w_ss); 
r   =log(r_ss); 
y   =log(y_ss); 
kap =log(kap_ss); 
innv=log(inv_ss); 
z   =log(z_ss);
g   =log(g_ss);
end;

shocks;

var e_z; stderr 0.01;
var e_g; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order = 1, hp_filter=1600, nograph);
