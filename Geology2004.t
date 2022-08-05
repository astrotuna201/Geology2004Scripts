
/*
_path="Geology2004_AnalysisScript/";
include Geology2004.t;
*/

macro main{
%setup_gnuplot;
%load_all_data;
%eccentricity_demodulation;
%obliquity_demodulation;
%precession_demodulation;
%prececc_demodulation;
};




macro setup_gnuplot{
/* set up gnuplot env.
_graph=grace;
*/
_graph=gnuplot;
opt="notitle w l";
gnuplot;
reset;
set term aqua;
set mxtics 5;
end;
};


macro load_all_data{
%make_etp_data;
%load_154_data;
%load_199_data;
};

macro load_etp_data{
_path="trip_analysis/";
file1="etp_la1993_adj.txt";
file2="etp_la2003_adj.txt";
tmin=16;
tmax=32;


vnumR	t93_init,etp93_init;	/* original data */
vnumR	t93_clip,etp93_clip;	/* clipped (selected time) data */
vnumR	t93_ip,etp93_ip;		/* interpolated data */

vnumR	t03_init,etp03_init;	/* original data */
vnumR	t03_clip,etp03_clip;	/* clipped (selected time) data */
vnumR	t03_ip,etp03_ip;		/* interpolated data */

vnumR	tbase,tbase_num;		/* used to i.p. to nearest base 2 */


read(file1,t93_init,etp93_init);
read(file2,t03_init,etp03_init);

/* clip to selected time interval (tmin,tmax global vars) */

t93_clip=select((t93_init <=tmax)&&(t93_init >=tmin),t93_init);
t03_clip=select((t03_init <=tmax)&&(t03_init >=tmin),t03_init);
etp93_clip=select((t93_init <= tmax)&&(t93_init >= tmin),etp93_init);
etp03_clip=select((t03_init <= tmax)&&(t03_init >= tmin),etp03_init);

/* interpolate to nearest base 2 */
tbase=log(size(t93_clip))/log(2);
tbase_num=2^(int(tbase)+1)-1;
t93_ip=min(t93_clip),max(t93_clip),(max(t93_clip)-min(t93_clip))/tbase_num;
etp93_ip=interpol(LINEAR,t93_clip,etp93_clip,t93_ip);

tbase=log(size(t03_clip))/log(2);
tbase_num=2^(int(tbase)+1)-1;
t03_ip=min(t03_clip),max(t03_clip),(max(t03_clip)-min(t03_clip))/tbase_num;
etp03_ip=interpol(LINEAR,t03_clip,etp03_clip,t03_ip);
};

macro make_etp_data{
_path="trip_analysis/";
file1="CLIVARN1993.asc";
file2="CLIVARN2003.asc";
vnumR	t93_ip,t03_ip,etp93_ip,etp03_ip;
vnumR	t93_clip,t03_clip,etp93_clip,etp03_clip;
vnumR	etp93,etp03;
vnumR	t93,e93,o93,p93;
vnumR	t03,e03,o03,p03;
vnumR	t_e_93,t_o_93,t_p_93;
vnumR	t_e_03,t_o_03,t_p_03;
vnumR	Eampl,Oampl,Pampl;
vnumR	Elag93,Elag03,Olag93,Olag03,Plag93,Plag03;
vnumR tmin,tmax;
vnumR	tbase,tbase_num;		/* used to i.p. to nearest base 2 */

tmin=17;
tmax=30;
Eampl=1.0;
Oampl=1.987549658;
Pampl=0.733508485;
Elag93=0;
Olag93=7.210971471/1e3; /*eq to 1.10507 angle if susref flipped (hi sus to low insolation)*/
Plag93=0;
Elag03=Elag93;
Olag03=Olag93;
Plag03=Plag93;

read(file1,t93,e93,o93,p93);
read(file2,t03,e03,o03,p03);
/*adjust times to Ma, 93 dat are pos, 03 are neg */
t93=t93/1000.0;
t03=t03/-1000.0;
/*apply time lags */
t_e_93=t93-Elag93;
t_o_93=t93-Olag93;
t_p_93=t93-Plag93;
t_e_03=t93-Elag03;
t_o_03=t93-Olag03;
t_p_03=t93-Plag03;
/* normalise all data */
e93=(e93-%avg[e93])/%stdev[e93];
o93=(o93-%avg[o93])/%stdev[o93];
p93=(p93-%avg[p93])/%stdev[p93];
e03=(e03-%avg[e03])/%stdev[e03];
o03=(o03-%avg[o03])/%stdev[o03];
p03=(p03-%avg[p03])/%stdev[p03];
/* re-interpolate lagged data */
e93=interpol(LINEAR,t_e_93,e93,t93);
o93=interpol(LINEAR,t_o_93,o93,t93);
p93=interpol(LINEAR,t_p_93,p93,t93);
e03=interpol(LINEAR,t_e_03,e03,t03);
o03=interpol(LINEAR,t_o_03,o03,t03);
p03=interpol(LINEAR,t_p_03,p03,t03);
/* create weighted ETP data */
etp93=e93*Eampl+o93*Oampl+p93*Pampl;
etp03=e03*Eampl+o03*Oampl+p03*Pampl;
etp93=(etp93-%avg[etp93])/%stdev[etp93];
etp03=(etp03-%avg[etp03])/%stdev[etp03];
/* clip to desired tmin, tmax range */
t93_clip=select((t93 <=tmax)&&(t93 >=tmin),t93);
t03_clip=select((t03 <=tmax)&&(t03 >=tmin),t03);
etp93_clip=select((t93 <= tmax)&&(t93 >= tmin),etp93);
etp03_clip=select((t03 <= tmax)&&(t03 >= tmin),etp03);
/* interpolate to nearest base 2 */
tbase=log(size(t93_clip))/log(2);
tbase_num=2^(int(tbase)+1)-1;
t93_ip=min(t93_clip),max(t93_clip),(max(t93_clip)-min(t93_clip))/tbase_num;
etp93_ip=interpol(LINEAR,t93_clip,etp93_clip,t93_ip);

tbase=log(size(t03_clip))/log(2);
tbase_num=2^(int(tbase)+1)-1;
t03_ip=min(t03_clip),max(t03_clip),(max(t03_clip)-min(t03_clip))/tbase_num;
etp03_ip=interpol(LINEAR,t03_clip,etp03_clip,t03_ip);
write(ETP1993adj_clip.txt,t93_ip,etp93_ip);
write(ETP2003adj_clip.txt,t03_ip,etp03_ip);
};

macro load_154_data{
_path="trip_analysis/";
file1="susref1993_sm9.txt";
file2="susref2003_smooth3b4.txt";
tmin=17;
tmax=30;


vnumR   t_sr93_init,sr93_init;	/* original data */
vnumR   t_sr93_clip,sr93_clip;	/* clipped (selected time) data */
vnumR   t_sr93_ip,sr93_ip;		/* interpolated data */
vnumR   sr93_ip_sm;				/* smoothed with 300pt mov avg*/

vnumR   t_sr03_init,sr03_init;	/* original data */
vnumR   t_sr03_clip,sr03_clip;	/* clipped (selected time) data */
vnumR   t_sr03_ip,sr03_ip;		/* interpolated data */
vnumR   sr03_ip_sm;				/* smoothed with 300pt mov avg*/
vnumR   tbase,tbase_num;		/* used to i.p. to nearest base 2 */


read(file1,t_sr93_init,sr93_init);
read(file2,t_sr03_init,sr03_init);

/* clip to selected time interval (tmin,tmax global vars) */

t_sr93_clip=select((t_sr93_init <=tmax)&&(t_sr93_init >=tmin),t_sr93_init);
t_sr03_clip=select((t_sr03_init <=tmax)&&(t_sr03_init >=tmin),t_sr03_init);
sr93_clip=select((t_sr93_init <= tmax)&&(t_sr93_init >= tmin),sr93_init);
sr03_clip=select((t_sr03_init <= tmax)&&(t_sr03_init >= tmin),sr03_init);

/* interpolate to nearest base 2 */
tbase=log(size(t_sr93_clip))/log(2);
tbase_num=2^(int(tbase)+1)-1;
t_sr93_ip=min(t_sr93_clip),max(t_sr93_clip),(max(t_sr93_clip)-min(t_sr93_clip))/tbase_num;
sr93_ip=interpol(LINEAR,t_sr93_clip,sr93_clip,t_sr93_ip);
/*sr93_ip_sm=sr93_ip-%smooth2[300,sr93_ip];*/
sr93_ip_sm=(sr93_ip-%avg[sr93_ip])/%stdev[sr93_ip];
/*sr93_ip_sm=sr93_ip;*/

tbase=log(size(t_sr03_clip))/log(2);
tbase_num=2^(int(tbase)+1)-1;
t_sr03_ip=min(t_sr03_clip),max(t_sr03_clip),(max(t_sr03_clip)-min(t_sr03_clip))/tbase_num;
sr03_ip=interpol(LINEAR,t_sr03_clip,sr03_clip,t_sr03_ip);
/*sr03_ip_sm=sr03_ip-%smooth2[300,sr03_ip];*/
sr03_ip_sm=(sr03_ip-%avg[sr03_ip])/%stdev[sr03_ip];
/*sr03_ip_sm=sr03_ip;*/
};


macro load_199_data{
_path="trip_analysis/";
file1="LA1218_depth.txt";
/*map_file="LA_good100ky_ptrs4.txt";
map_file="LA100ky_ptr6.txt";
map_file="depth_age_njsck.txt";
map_file="LA100ky_ptr9.txt";*/
map_file="LA100ky_ptr11.txt";
tmin=21.0;
tmax=28.6;
vnumR	dep_199,dat_199;			/* original data vs. depth */
vnumR	t_199_init,dat_199_init;	/* original data */
vnumR	t_199_clip,dat_199_clip;	/* clipped (selected time) data */
vnumR	t_190_ip,dat_199_ip;		/* interpolated data */

/* read in data vs. depth */
read(file1,dep_199,dat_199);
/* read in depth to time mapping function */
vnumR z_base, t_base; /* mapping depths and times*/
read(map_file,z_base,t_base);
/* interpolate depth to time */
t_199_init=interpol(LINEAR,z_base,t_base,dep_199);
dat_199_init=dat_199;
/* clip to selected time interval (tmin,tmax global vars) */
t_199_clip=select((t_199_init <=tmax)&&(t_199_init >=tmin),t_199_init);
dat_199_clip=select((t_199_init <= tmax)&&(t_199_init >= tmin),dat_199_init);
/* interpolate to nearest base 2 */
tbase=log(size(t_199_clip))/log(2);
tbase_num=2^(int(tbase)+1)-1;
t_199_ip=min(t_199_clip),max(t_199_clip),(max(t_199_clip)-min(t_199_clip))/tbase_num;
dat_199_ip=interpol(LINEAR,t_199_clip,dat_199_clip,t_199_ip);
dat_199_ip_sm=dat_199_ip-%smooth2[600,dat_199_ip];
dat_199_ip_sm=(dat_199_ip_sm-%avg[dat_199_ip_sm])/%stdev[dat_199_ip_sm];
};

macro demod[demod_freq,demod_t_in,demod_dat_in,fmin,fmax]{

vnumR demodul,demodulabs,yamp,yphase;
dtour=2*pi;
demodul=demod_dat_in*exp(-i*dtour*demod_freq*demod_t_in);
demodabs=real(demodul);
%filtreC[demodul,demod_t_in,1,fmin,fmax];
yphase=atan2(_fi_x,_fi_y);
yamp=2*abs(_fi_x+i*_fi_y);
return(yamp);
};


macro demodphase[demod_freq,demod_t_in,demod_dat_in,fmin,fmax]{

vnumR demodul,demodulabs,yamp,yphase;
dtour=2*pi;
demodul=demod_dat_in*exp(-i*dtour*demod_freq*demod_t_in);
demodabs=real(demodul);
%filtreC[demodul,demod_t_in,1,fmin,fmax];
yphase=atan2(_fi_x,_fi_y);
yamp=2*abs(_fi_x+i*_fi_y);
return(yphase);
};




macro precession_demodulation{
vnumR etp93_pdem,etp03_pdem;
vnumR sr93_pdem,sr03_pdem;
vnumR sr93prec,sr93prec_hil,sr93_pdem_frac,sr93_pdem_frac_cut;
vnumR sr03prec,sr03prec_hil,sr03_pdem_frac,sr03_pdem_frac_cut;
vnumR etp93prec,etp93prec_hil,etp93_pdem_frac;
vnumR etp03prec,etp03prec_hil,etp03_pdem_frac;
vnumR sr93_pdem_frac_cut,sr03_pdem_frac_cut,t_sr93_ip_cut,t_sr03_ip_cut;
vnumR dem_lopass;

/* demodulate ETP target and Leg 154 data at frequencies
identified with subnafr */
/*sr93_ip_sm:
"/y	years	Ampl	phase	cyc/My
68.378662	18953	0.038957	24.18	52.76209571
68.9037	18809	0.047049	-7.993	53.16603754
*/
/*sr03_ip_sm:
"/y	years	Ampl	phase	cyc/My
68.472626	18927	0.06839	56.346	52.83457495
69.007755	18780	0.066377	-42.518	53.24813632
*/
/*etp03_ip:
"/y	years	Ampl	phase	cyc/My
68.479607	18925	0.067058	13.435	52.84015852
69.019163	18777	0.095688	-121.403	53.25664377
*/

dem_lopass=0.9;
etp93_pdem=((%demod[53.16038488,t93_ip,etp93_ip,-dem_lopass,dem_lopass]+%demod[52.71626389,t93_ip,etp93_ip,-dem_lopass,dem_lopass]))/2;
etp03_pdem=((%demod[53.25664377,t03_ip,etp03_ip,-dem_lopass,dem_lopass]+%demod[52.84015852,t03_ip,etp03_ip,-dem_lopass,dem_lopass]))/2;
/*
sr93_pdem=((%demod[53.16038488,t_sr93_ip,sr93_ip_sm,-dem_lopass,dem_lopass]+%demod[52.71626389,t_sr93_ip,sr93_ip_sm,-dem_lopass,dem_lopass]))/2;
sr03_pdem=((%demod[53.25664377,t_sr03_ip,sr03_ip_sm,-dem_lopass,dem_lopass]+%demod[52.84015852,t_sr03_ip,sr03_ip_sm,-dem_lopass,dem_lopass]))/2;
*/
sr93_pdem=((%demod[53.16603754,t_sr93_ip,sr93_ip_sm,-dem_lopass,dem_lopass]+%demod[52.76209571,t_sr93_ip,sr93_ip_sm,-dem_lopass,dem_lopass]))/2;
sr03_pdem=((%demod[53.24813632,t_sr03_ip,sr03_ip_sm,-dem_lopass,dem_lopass]+%demod[52.83457495,t_sr03_ip,sr03_ip_sm,-dem_lopass,dem_lopass]))/2;


/*
1/53.16038488;
18811 ky
1/52.71626389;
18969 ky
1/53.25664377;
18776 ky
1/52.84015852;
18925 ky
*/
/* 	climatic precession signal in 154 data varies in total amplitude,
	hence calculate amplitude modulation of ~19ky cycle as fraction
	of total precession amplitude */
	
%filtreR[sr93_ip_sm,t_sr93_ip,1,40,55];
sr93prec=_fi_x;
sr93prec_hil=abs(%hilbert[sr93prec]);
sr93prec_hil=%smooth2[2000,sr93prec_hil];
sr93_pdem_frac=sr93_pdem/sr93prec_hil;

%filtreR[sr03_ip_sm,t_sr03_ip,1,40,55];
sr03prec=_fi_x;
sr03prec_hil=abs(%hilbert[sr03prec]);
sr03prec_hil=%smooth2[2000,sr03prec_hil];
sr03_pdem_frac=sr03_pdem/sr03prec_hil;

%filtreR[etp93_ip,t93_ip,1,40,55];
etp93prec=_fi_x;
etp93prec_hil=abs(%hilbert[etp93prec]);
etp93prec_hil=%smooth2[2000,etp93prec_hil];
etp93_pdem_frac=etp93_pdem/etp93prec_hil;

%filtreR[etp03_ip,t03_ip,1,40,55];
etp03prec=_fi_x;
etp03prec_hil=abs(%hilbert[etp03prec]);
etp03prec_hil=%smooth2[2000,etp03prec_hil];
etp03_pdem_frac=etp03_pdem/etp03prec_hil;

/* smoothing and data processing affected geol data, 
which are only loaded to 30My - cut off to 29.5Ma, as
data quality is poorer there anyway */

tmin=17.5;
tmax=29.5;
t_sr93_ip_cut=select((t_sr93_ip <=tmax)&&(t_sr93_ip >=tmin),t_sr93_ip);
t_sr03_ip_cut=select((t_sr03_ip <=tmax)&&(t_sr03_ip >=tmin),t_sr03_ip);
sr93_pdem_frac_cut=select((t_sr93_ip<=tmax)&&(t_sr93_ip>=tmin),sr93_pdem_frac);
sr03_pdem_frac_cut=select((t_sr03_ip <= tmax)&&(t_sr03_ip >= tmin),sr03_pdem_frac);

plot(t_sr03_ip_cut,sr03_pdem_frac_cut+0.5,opt);
replot(t03_ip,etp03_pdem_frac+0.5,opt);
replot(t_sr93_ip_cut,sr93_pdem_frac_cut,opt);
replot(t93_ip,etp93_pdem_frac,opt);
write("astr_prec_93_03_res.txt",t03_ip,etp93_pdem_frac,etp03_pdem_frac);
write("geol_prec_93_03_res.txt",t_sr03_ip_cut,sr93_pdem_frac_cut,sr03_pdem_frac_cut);
write("astr_prec_93_03_pdem.txt",t03_ip,etp93_pdem,etp03_pdem);
write("geol_prec_93_03_pdem.txt",t_sr03_ip,sr93_pdem,sr03_pdem);
write("smoothed_154_93_03_data.txt",t_sr03_ip,sr93_ip_sm,sr03_ip_sm);
write("astr_data.txt",t03_ip,etp93_ip,etp03_ip);

};

macro obliquity_demodulation{
vnumR etp93_odem,etp03_odem;
vnumR sr93_odem,sr03_odem;
vnumR sr93obl,sr93obl_hil,sr93_odem_frac,sr93_odem_frac_cut;
vnumR sr03obl,sr03obl_hil,sr03_odem_frac,sr03_odem_frac_cut;
vnumR etp93obl,etp93obl_hil,etp93_odem_frac;
vnumR etp03obl,etp03obl_hil,etp03_odem_frac;
vnumR sr93_odem_frac_cut,sr03_odem_frac_cut,t_sr93_ip_cut,t_sr03_ip_cut;
vnumR sr03_odem_freq,sr93_odem_freq,etp03_odem_freq,etp93_odem_freq;
vnumR dem_lopass;
sr03_odem_freq=32.250643; /*from subnafr of sr03_ip_sm */
sr93_odem_freq=32.175073;
etp03_odem_freq=32.263647;
etp93_odem_freq=32.18497;

sr03_odem_freq=%asec_to_Ma[sr03_odem_freq];/* convert from arcsec/yr to cyc/My */
sr93_odem_freq=%asec_to_Ma[sr93_odem_freq];
etp03_odem_freq=%asec_to_Ma[etp03_odem_freq];
etp93_odem_freq=%asec_to_Ma[etp93_odem_freq];
/* frequencies in cycles/My
sr03_odem_freq = 24.8847554012346
sr93_odem_freq = 24.8264452160494
etp03_odem_freq = 24.8947893518519
etp93_odem_freq = 24.8340817901235
*/

dem_lopass=1.0;
etp03_odem=(%demod[etp03_odem_freq,t03_ip,etp03_ip,-dem_lopass,dem_lopass]);
etp93_odem=(%demod[etp93_odem_freq, t93_ip,etp93_ip,-dem_lopass,dem_lopass]);
sr03_odem=(%demod[sr03_odem_freq,t_sr03_ip,sr03_ip_sm,-dem_lopass,dem_lopass]);
sr93_odem=(%demod[etp93_odem_freq,t_sr93_ip,sr93_ip_sm,-dem_lopass,dem_lopass]);


etp03_odem_phase=(%demodphase[etp03_odem_freq,t03_ip,etp03_ip,-dem_lopass,dem_lopass]);
etp93_odem_phase=(%demodphase[etp93_odem_freq, t93_ip,etp93_ip,-dem_lopass,dem_lopass]);
sr03_odem_phase=(%demodphase[sr03_odem_freq,t_sr03_ip,sr03_ip_sm,-dem_lopass,dem_lopass]);
sr93_odem_phase=(%demodphase[etp93_odem_freq,t_sr93_ip,sr93_ip_sm,-dem_lopass,dem_lopass]);



/*	calculate amplitude modulation of ~19ky cycle as fraction
	of total precession amplitude */
	
%filtreR[sr93_ip_sm,t_sr93_ip,1,20,30];
sr93obl=_fi_x;
sr93obl_hil=abs(%hilbert[sr93obl]);
sr93obl_hil=%smooth2[2000,sr93obl_hil]*1.8/1.1;
sr93_odem_frac=sr93_odem/sr93obl_hil;

%filtreR[sr03_ip_sm,t_sr03_ip,1,20,30];
sr03obl=_fi_x;
sr03obl_hil=abs(%hilbert[sr03obl]);
sr03obl_hil=%smooth2[2000,sr03obl_hil]*1.8/1.1;
sr03_odem_frac=sr03_odem/sr03obl_hil;

%filtreR[etp93_ip,t93_ip,1,20,30];
etp93obl=_fi_x;
etp93obl_hil=abs(%hilbert[etp93obl]);
etp93obl_hil=%smooth2[2000,etp93obl_hil]*1.8/1.1;
etp93_odem_frac=etp93_odem/etp93obl_hil;

%filtreR[etp03_ip,t03_ip,1,20,30];
etp03obl=_fi_x;
etp03obl_hil=abs(%hilbert[etp03obl]);
etp03obl_hil=%smooth2[2000,etp03obl_hil]*1.8/1.1;
etp03_odem_frac=etp03_odem/etp03obl_hil;

tmin=17.5;
tmax=30;
t_sr93_ip_cut=select((t_sr93_ip <=tmax)&&(t_sr93_ip >=tmin),t_sr93_ip);
t_sr03_ip_cut=select((t_sr03_ip <=tmax)&&(t_sr03_ip >=tmin),t_sr03_ip);
sr93_odem_frac_cut=select((t_sr93_ip<=tmax)&&(t_sr93_ip>=tmin),sr93_odem_frac);
sr03_odem_frac_cut=select((t_sr03_ip<=tmax)&&(t_sr03_ip>=tmin),sr03_odem_frac);

plot(t_sr03_ip_cut,sr03_odem_frac_cut+1.0,opt);
replot(t_sr93_ip_cut,sr93_odem_frac_cut,opt);
replot(t03_ip,etp03_odem_frac+1.0,opt);
replot(t93_ip,etp93_odem_frac,opt);
write("astr_obl_93_03_res.txt",t03_ip,etp93_odem_frac,etp03_odem_frac);
write("geol_obl_93_03_res.txt",t_sr03_ip_cut,sr93_odem_frac_cut,sr03_odem_frac_cut);
write("astr_obl_93_03_odem.txt",t03_ip,etp93_odem,etp03_odem);
write("geol_obl_93_03_odem.txt",t_sr03_ip,sr93_odem,sr03_odem);

};


macro prececc_demodulation{
/* demodulate ETP target and Leg 154 data at strongest precession frequency
identified with subnafr, and then demodulate again for short and long eccentricity */
vnumR etp03_pdem,etp93_pdem,sr03_pdem,sr93_pdem;/* precession demodulations */
vnumR etp03prec,etp93prec,sr03prec,sr93prec;/* filter around precession band */
vnumR etp03prec_hil,etp93prec_hil,sr03prec_hil,sr93prec_hil;/* ampl of filter*/
vnumR etp03_pdem_frac,etp93_pdem_frac,sr03_pdem_frac,sr93_pdem_frac;/* precession demodulations normalised */
vnumR w1,w2,w3,w4,w5,wsum; /* weights used for ind. demods of prec */
vnumR w1dem,w2dem,w3dem,w4dem,w5dem, w23dem,w45dem;
vnumR f1,f2,f3,f4,f5,f23,f45; /* frequencies for ecc at which to demod */
vnumR fmed03,fmed93,fmedsr03,fmedsr93;
vnumR f1_03,f1_93,f23_03,f23_93,f45_03,f45_93;
vnumR final_filter_freq;
vnumR precdemfreq,filtwidth,filtwidth2;

final_filter_freq=1.2;
filtwidth2=0.8;

/* indices:
1	406ky
2	94ky
3	98ky
4	124ky
5	131ky
*/
/****************DEMDULATION FOR ETP2003 *******************/
/* etp03 strongest precession frequency */
NTERM=10;
/* strongest prec freq det by subnafr 
55.368147	23407	0.145492	-90.864		42.72226257
58.562675	22130	0.107499	118.43		45.18752824
68.479607	18925	0.067058	13.435		52.84015852
69.019163	18777	0.095688	-121.403	53.25664377
*/
/*fmed03=(55.368147*0.145492+58.562675*0.107499+68.479607*0.067058+69.019163*0.095688)/
(0.145492+0.107499+0.067058+0.095688);*/

fmed03=55.368147;
fmed03=%asec_to_Ma[fmed03];
/*fmed03 = 61.4510197516675 "/a *//* ca 21.09 ky */
/*fmed03 = 47.4159103022126 1/My*/
precdemfreq=fmed03;

filtwidth=12;
%filtreR[etp03_ip,t03_ip,1,40,58];
etp03prec=_fi_x;
etp03_pdem=%demod[precdemfreq,t03_ip,etp03prec,-filtwidth,filtwidth];
/* now normalise to strength of precession signal */
etp03prec_hil=abs(%hilbert[etp03prec]);
etp03prec_hil=%smooth2[2000,etp03prec_hil];
etp03_pdem_frac=etp03_pdem/etp03prec_hil;
etp03_pdem=(etp03_pdem_frac-%avg[etp03_pdem_frac])/%stdev[etp03_pdem_frac];

/*plot(t03_ip,etp03_pdem,opt);*/
/*%subnafr[t03_ip*1e6,etp03_pdem];*/
/*
1        3.195583        405560      0.388599   -159.449  
2       -3.195583       -405560      0.388599    159.451  
3      -13.652271        -94929      0.311791     44.569  
4       13.652271         94929      0.311791    -44.569  
5       10.458025        123924      0.243787    112.717  
6      -10.458025       -123924      0.243787   -112.717  
7       13.108706         98866      0.212813    115.017  
8      -13.108706        -98866      0.212813   -115.017  
9       -9.915746       -130701      0.161474    103.627  
10        9.915746        130701      0.161474   -103.627  
*/
/*
f1=tf[1];f2=tf[4];f3=tf[7];f4=tf[5];f5=tf[10];w1=abs(za[1]);w2=abs(za[4]);w3=abs(za[7]);w4=abs(za[5]);w5=abs(za[10]);*/
/* these are in arcsec/a - need to convert to Ma */
f1 = 3.19558254199462;
f2 = 13.6522705586679;
f3 = 13.1087061126672;
f4 = 10.4580254020956;
f5 = 9.91574570532256;
w1 = 0.388598607459231;
w2 = 0.31179104791276;
w3 = 0.21281312734302;
w4 = 0.243786858319462;
w5 = 0.161473561826929;
f1=%asec_to_Ma[f1];f2=%asec_to_Ma[f2];f3=%asec_to_Ma[f3];f4=%asec_to_Ma[f4];f5=%asec_to_Ma[f5];


wsum=w1+w2+w3+w4+w5;
/*filtwidth2=0.9;*/

f1_03=f1;
f23_03=(f2+f3)/2;
f45_03=(f4+f5)/2;
/*
f1 = 2.46566111764145;
f23 = 10.324487658334;
f45 = 7.86000502810213;
*/

w1dem=w1/wsum*%demod[f1,t03_ip,etp03_pdem,-filtwidth2,filtwidth2];
w1dem=2*%avg[w1dem]-w1dem; /* flip 406 ky component */
w23dem=(w2+w3)/wsum*%demod[f23_03,t03_ip,etp03_pdem,-filtwidth2,filtwidth2];
w45dem=(w4+w5)/wsum*%demod[f45_03,t03_ip,etp03_pdem,-filtwidth2,filtwidth2];

w2dem=w2/wsum*%demod[f2,t03_ip,etp03_pdem,-filtwidth2,filtwidth2];
w3dem=w3/wsum*%demod[f3,t03_ip,etp03_pdem,-filtwidth2,filtwidth2];
w4dem=w4/wsum*%demod[f4,t03_ip,etp03_pdem,-filtwidth2,filtwidth2];
w5dem=w5/wsum*%demod[f5,t03_ip,etp03_pdem,-filtwidth2,filtwidth2];
etp03_prececc_dem=w1dem+w2dem+w3dem+w4dem+w5dem;
etp03_prececc_dem=w1dem+w23dem+w45dem;

%filtreR[etp03_prececc_dem,t03_ip,1,-final_filter_freq,final_filter_freq];
etp03_prececc_dem=_fi_x;
/*
plot(t03_ip,etp03_prececc_dem,opt);
replot(t03_ip,w1dem,opt);
replot(t03_ip,w23dem,opt);replot(t03_ip,w45dem,opt);

replot(t03_ip,w2dem,opt);replot(t03_ip,w3dem,opt);
replot(t03_ip,w4dem,opt);replot(t03_ip,w5dem,opt);*/

/************END DEMDULATION FOR ETP2003 *******************/


/****************DEMDULATION FOR ETP1993 *******************/

/* etp93 strongest precession frequency */
NTERM=10;
/* strongest prec freq det by subnafr 
55.287686	23441	0.152767	146.828	42.66029606
58.495727	22155	0.123971	-22.624	45.13653803
68.245266	18990	0.055411	100.273	52.65929437
68.395733	18949	0.055894	-81.762	52.77323342
68.895145	18811	0.117171	29.958	53.16038488
*/
fmed93=(55.287686*0.152767+58.495727*0.123971+68.245266*0.055411+68.395733*0.055894+68.895145*0.117171)/
(0.152767+0.123971+0.055411+0.055894+0.117171);
fmed93=55.287686;

fmed93=%asec_to_Ma[fmed93];

/*fmed93 = 62.1021394270586 "/a *//* ca 21.09 ky */
/*fmed93 = 47.9183174591502 1/My*/
precdemfreq=fmed93;

filtwidth=12;
%filtreR[etp93_ip,t93_ip,1,40,58];
etp93prec=_fi_x;
etp93_pdem=%demod[precdemfreq,t93_ip,etp93prec,-filtwidth,filtwidth];
/*plot(t93_ip,etp93_pdem,opt);*/
/* now normalise to strength of precession signal */
etp93prec_hil=abs(%hilbert[etp93prec]);
etp93prec_hil=%smooth2[2000,etp93prec_hil];
etp93_pdem_frac=etp93_pdem/etp93prec_hil;

etp93_pdem=(etp93_pdem_frac-%avg[etp93_pdem_frac])/%stdev[etp93_pdem_frac];
/*plot(t93_ip,etp93_pdem,opt);*/
/*%subnafr[t93_ip*1e6,etp93_pdem];*/
/*
1       -3.207530       -404049      0.398342    168.298  
2        3.207529        404049      0.398342   -168.297  
3       13.090340         99004      0.166584   -131.970  
4      -13.090340        -99004      0.166584    131.970  
5       -9.880871       -131163      0.145161    -48.468  
6        9.880871        131163      0.145161     48.468  
7      -13.619869        -95155      0.302092   -157.930  
8       13.619869         95155      0.302092    157.930  
9      -10.413777       -124451      0.253614     45.433  
10       10.413777        124451      0.253614    -45.433  
       */
/*f1=tf[2];f2=tf[8];f3=tf[3];f4=tf[10];f5=tf[6];w1=abs(za[2]);w2=abs(za[8]);w3=abs(za[3]);w4=abs(za[10]);w5=abs(za[6]);*/
 /* these are in arcsec/a - need to convert to Ma */
f1 = 3.20752947655787;
f2 = 13.6198694742647;
f3 = 13.0903403499093;
f4 = 10.4137773782909;
f5 = 9.8808705526489;
w1 = 0.398341718960711;
w2 = 0.302091997266692;
w3 = 0.166584199905735;
w4 = 0.253614448769277;
w5 = 0.145161158037453;
f1=%asec_to_Ma[f1];f2=%asec_to_Ma[f2];f3=%asec_to_Ma[f3];f4=%asec_to_Ma[f4];f5=%asec_to_Ma[f5];

wsum=w1+w2+w3+w4+w5;
/*filtwidth2=0.9;*/

f1_93=f1;
f23_93=(f2+f3)/2;
f45_93=(f4+f5)/2;
/*
f1 = 2.47497166947414;
f23 = 10.3049122463161;
f45 = 7.82958153157725;
*/

w1dem=w1/wsum*%demod[f1,t93_ip,etp93_pdem,-filtwidth2,filtwidth2];
w1dem=2*%avg[w1dem]-w1dem; /* flip 406 ky component */
w23dem=(w2+w3)/wsum*%demod[f23_93,t93_ip,etp93_pdem,-filtwidth2,filtwidth2];
w45dem=(w4+w5)/wsum*%demod[f45_93,t93_ip,etp93_pdem,-filtwidth2,filtwidth2];

w2dem=w2/wsum*%demod[f2,t93_ip,etp93_pdem,-filtwidth2,filtwidth2];
w3dem=w3/wsum*%demod[f3,t93_ip,etp93_pdem,-filtwidth2,filtwidth2];
w4dem=w4/wsum*%demod[f4,t93_ip,etp93_pdem,-filtwidth2,filtwidth2];
w5dem=w5/wsum*%demod[f5,t93_ip,etp93_pdem,-filtwidth2,filtwidth2];

etp93_prececc_dem=w1dem+w2dem+w3dem+w4dem+w5dem;
etp93_prececc_dem=w1dem+w23dem+w45dem;

%filtreR[etp93_prececc_dem,t93_ip,1,-final_filter_freq,final_filter_freq];
etp93_prececc_dem=_fi_x;



/*plot(t93_ip,w1dem,opt);
replot(t93_ip,w2dem,opt);replot(t93_ip,w3dem,opt);
replot(t93_ip,w4dem,opt);replot(t93_ip,w5dem,opt);*/
/************END DEMDULATION FOR ETP1993 *******************/
/*plot(t03_ip,etp03_prececc_dem,"title 'etp03_prececc' w l");
replot(t93_ip,etp93_prececc_dem,"title 'etp93_prececc' w l");*/
write("astr_prececc_93_03_res.txt",t03_ip,etp93_prececc_dem,etp03_prececc_dem);


/****************DEMDULATION FOR SR2003 *******************/
/* sr03 strongest precession frequency */
NTERM=100;
precdemfreq=fmed03; /* determined with subnafr*/
/* actual values for strongest prec freq of sr03
	55.345996	23416	0.129053	37.03	42.70584216
	58.535232	22141	0.097262	-71.174	45.16507836
	68.472626	18927	0.06839	56.346	52.83457495
	69.007755	18780	0.066377	-42.518	53.24813632
*/
/*fmedsr03=(55.340951*0.109128+58.542815*0.07831+68.46177*0.052892+68.995968*0.054302)/
(0.109128+0.07831+0.052892+0.054302);*/

fmedsr03=55.345996;
fmedsr03=%asec_to_Ma[fmedsr03];
/*fmedsr03 = 61.0640871886082 "/a *//* ca 21.09 ky */
/*fmedsr03 = 47.1173512257779 1/My*/
precdemfreq=fmedsr03; /* determined with subnafr*/
filtwidth=12;
%filtreR[sr03_ip_sm,t_sr03_ip,1,40,58];
sr03prec=_fi_x;
sr03_pdem=%demod[precdemfreq,t_sr03_ip,sr03prec,-filtwidth,filtwidth];
/*%filtreR[sr03_pdem,t_sr03_ip,1,1,12];
sr03_pdem=_fi_x;*/
/* now normalise to strength of precession signal */
sr03prec_hil=abs(%hilbert[sr03prec]);
sr03prec_hil=%smooth2[2000,sr03prec_hil];
sr03_pdem_frac=sr03_pdem/sr03prec_hil;

sr03_pdem=(sr03_pdem_frac-%avg[sr03_pdem_frac])/%stdev[sr03_pdem_frac];

/*plot(t03_ip,etp03_pdem,"title 'etp03_pdem' w l");
replot(t_sr03_ip,sr03_pdem,"title 'sr03_pdem' w l");*/

/*%subnafr[t_sr03_ip*1e6,sr03_pdem];*/
/*
1      -13.658043        -94889      0.203399     60.920  
2       13.658043         94889      0.203399    -60.920  
3       -3.193258       -405855      0.320273    144.205  
4        3.193258        405855      0.320272   -144.207  
5       13.127648         98723      0.175188     18.952  
6      -13.127648        -98723      0.175188    -18.952  
7      -10.459632       -123905      0.142474   -115.067  
8       10.459632        123905      0.142474    115.067  
9       -6.427815       -201624      0.116330    -71.617  
10        6.427815        201624      0.116330     71.617  
11       -1.420615       -912281      0.162871   -108.361  
12        1.420615        912281      0.162873    108.359  
13        2.612199        496134      0.122813    133.915  
14       -2.612198       -496134      0.122810   -133.917  
15        3.889070        333242      0.068428     21.440  
16       -3.889070       -333242      0.068427    -21.440  
17        0.504543       2568659      0.092385    119.137  
18        0.721401       1796505      0.069245    -35.439  
19       -0.504572      -2568512      0.096627   -114.285  
20       -0.721307      -1796740      0.069512     34.485  
21       -5.636802       -229918      0.063229   -140.942  
22        5.636802        229918      0.063229    140.942  
23       11.692391        110841      0.065146     23.413  
24      -11.692391       -110841      0.065146    -23.413  
25        7.062361        183508      0.108561   -108.903  
26       -7.062361       -183508      0.108560    108.901  
27        3.474149        373041      0.111227     96.799  
28       -3.474148       -373041      0.111227    -96.801  
29        3.638244        356216      0.071678     68.792  
30       -3.638244       -356216      0.071678    -68.794  
31      -13.477922        -96157      0.057855     67.102  
32       -9.203343       -140818      0.065357    124.901  
33       13.477922         96157      0.057855    -67.102  
34        9.203343        140818      0.065357   -124.901  
35       -2.036373       -636426      0.076272   -117.578  
36        2.036374        636425      0.076269    117.567  
37       -0.171771      -7544950      0.043441   -167.307  
38        2.882887        449549      0.047059    146.557  
39       -2.882880       -449550      0.047060   -146.602  
40        0.350692       3695550      0.062520   -173.077  
41        1.096061       1182416      0.072245    178.523  
42       -5.422045       -239024      0.094607    -48.074  
43        5.422046        239024      0.094607     48.072  
44       -1.096188      -1182279      0.072189   -177.711  
45       16.357148         79231      0.056072   -140.784  
46      -16.357148        -79231      0.056072    140.784  
47        5.149209        251689      0.045626    103.375  
48       -5.149211       -251689      0.045626   -103.366  
49       -7.982243       -162360      0.065857    -12.256  
50        7.982243        162360      0.065857     12.257  
51       -1.587785       -816232      0.044085    -36.454  
52        1.587801        816223      0.044093     36.356  
53       12.830009        101013      0.047502   -106.874  
54      -12.830008       -101013      0.047502    106.873  
55        9.911752        130754      0.116609    -56.275  
56        9.751238        132906      0.073679    -34.482  
57       -7.300640       -177519      0.046068    142.281  
58       -9.911750       -130754      0.116608     56.262  
59       -9.751238       -132906      0.073678     34.480  
60       -0.346433      -3740982      0.063214    142.972  
61        8.525474        152015      0.040931    134.720  
62        7.300633        177519      0.046076   -142.234  
63       -8.525489       -152015      0.040943   -134.622  
64      -16.836340        -76976      0.045434   -109.013  
65       16.836340         76976      0.045434    109.013  
66        4.504358        287721      0.058904    142.451  
67       -4.504358       -287721      0.058904   -142.452  

*/
/*f1=tf[4];f2=tf[2];f3=tf[5];f4=tf[8];f5=tf[55];w1=abs(za[4]);w2=abs(za[2]);w3=abs(za[5]);w4=abs(za[8]);w5=abs(za[55]);*/
/*%subnafres[t_sr03_ip*1e6];*/
 /*
     -10.034461       -129155      0.040875    -97.668  
       9.942837        130345      0.130449     73.467  
       9.690782        133735      0.094723    -50.061  
 f5=tf[2];w5=abs(za[2]);     */
/* these are in arcsec/a - need to convert to Ma */

/* Frequencies from data*/ 
f1 = 3.19325840785723;
f2 = 13.6580434123333;
f3 = 13.1276483294571;
f4 = 10.4596315906993;
f5 = 9.9117518541922;
w1 = 0.320272059778701;
w2 = 0.203398820865281;
w3 = 0.175188433409666;
w4 = 0.142474159261663;
w5 = 0.116608515141883;

/* frequencies from ETP2003 */
/*f1 = 3.20087910542795;
f2 = 13.6693285057481;
f3 = 13.1176584138709;
f4 = 10.4519663458657;
f5 = 9.94533519422196;
w1 = 0.27434803586096;
w2 = 0.178534014211382;
w3 = 0.155580340260222;
w4 = 0.138708822704022;
w5 = 0.106367335885613;*/

f1=%asec_to_Ma[f1];f2=%asec_to_Ma[f2];f3=%asec_to_Ma[f3];f4=%asec_to_Ma[f4];f5=%asec_to_Ma[f5];

wsum=w1+w2+w3+w4+w5;
/*filtwidth2=0.9;*/

w1dem=w1/wsum*%demod[f1,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];
w1dem=2*%avg[w1dem]-w1dem; /* flip 406 ky component */
w2dem=w2/wsum*%demod[f2,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];
w3dem=w3/wsum*%demod[f3,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];
w4dem=w4/wsum*%demod[f4,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];
w5dem=w5/wsum*%demod[f5,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];
sr03_prececc_dem=w1dem+w2dem+w3dem+w4dem+w5dem;
/*%filtreR[sr03_prececc_dem,t_sr03_ip,1,0.3,0.7];
sr03_prececc_dem=_fi_x;*/
/*filtwidth2=0.8;*/
filtwidth3=0.8;
f23=(f2+f3)/2;f45=(f4+f5)/2;
w1dem=w1/wsum*%demod[f1,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];
w1dem=2*%avg[w1dem]-w1dem; /* flip 406 ky component */
w23dem=(w2+w3)/wsum*%demod[f23,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];
w45dem=(w4+w5)/wsum*%demod[f45,t_sr03_ip,sr03_pdem,-filtwidth2,filtwidth2];

sr03_prececc_dem=w1dem+w2dem+w3dem+w4dem+w5dem;
sr03_prececc_dem=w1dem+w23dem+w45dem;

%filtreR[sr03_prececc_dem,t_sr03_ip,1,-final_filter_freq,final_filter_freq];
sr03_prececc_dem=_fi_x;

/*
plot(t_sr03_ip,sr03_pdem/3,opt);
replot(t_sr03_ip,w1dem,opt);
replot(t_sr03_ip,w23dem,opt);
replot(t_sr03_ip,w45dem,opt);

plot(t_sr03_ip,w1dem,opt);replot(t_sr03_ip,w2dem,opt);replot(t_sr03_ip,w3dem,opt);
replot(t03_ip,w4dem,opt);replot(t_sr03_ip,w5dem,opt);*/
/*plot(t03_ip,etp03_prececc_dem,"title 'etp03_prececc_dem' w l");
replot(t93_ip,etp93_prececc_dem,"title 'etp93_prececc_dem' w l");
replot(t_sr03_ip,sr03_prececc_dem,"title 'sr03_prececc_dem' w l");*/

/************END DEMDULATION FOR SR2003 *******************/

/****************DEMDULATION FOR SR1993 *******************/
/* sr93 strongest precession frequency */
NTERM=100;
precdemfreq=fmed93; /* determined with subnafr*/
/* actual values for strongest prec freq of sr93
55.276038	23446	0.100563	-146.668	42.6511985
58.492445	22157	0.089724	-8.477	45.13246378
68.378662	18953	0.038957	24.18	52.76209571
68.9037	18809	0.047049	-7.993	53.16603754
*/
/*fmedsr93=(55.278883*0.091641+58.483049*0.057136+68.39095*0.037103+68.895004*0.03438)/
(0.091641+0.057136+0.037103+0.03438);*/
fmedsr93=55.276038;
fmedsr93=%asec_to_Ma[fmedsr93];
/*fmedsr03 = 60.4441081450876 "/a *//* ca 21.09 ky */
/*fmedsr03 = 46.6389723341726 1/My*/
precdemfreq=fmedsr93; /* determined with subnafr*/
fmedsr93=fmed93;
filtwidth=12;
%filtreR[sr93_ip_sm,t_sr93_ip,1,40,58];
sr93prec=_fi_x;
sr93_pdem=%demod[precdemfreq,t_sr93_ip,sr93prec,-filtwidth,filtwidth];
/*%filtreR[sr93_pdem,t_sr93_ip,1,0.8,12];
sr93_pdem=_fi_x;*/
/* now normalise to strength of precession signal */
sr93prec_hil=abs(%hilbert[sr93prec]);
sr93prec_hil=%smooth2[2000,sr93prec_hil];
sr93_pdem_frac=sr93_pdem/sr93prec_hil;

sr93_pdem=(sr93_pdem_frac-%avg[sr93_pdem_frac])/%stdev[sr93_pdem_frac];

/*
plot(t93_ip,etp93_pdem,"title 'etp93_pdem' w l");
replot(t_sr93_ip,sr93_pdem,"title 'sr93_pdem' w l");
*/

/*%subnafr[t_sr93_ip*1e6,sr93_pdem];*/
/*
1       -3.226656       -401654      0.349414    -77.927  
2        3.226656        401654      0.349414     77.926  
3      -10.079485       -128578      0.105049     62.469  
4       10.079485        128578      0.105048    -62.472  
5        0.529743       2446472      0.141586    -26.321  
6      -13.648037        -94959      0.121665    -29.755  
7       13.648037         94959      0.124520     15.739  
8       -0.529722      -2446568      0.141639     26.181  
9       -9.146517       -141693      0.103466    -87.310  
10      -10.415498       -124430      0.131429     43.720  
11       10.415479        124430      0.131427    -43.602  
12        9.146478        141694      0.103467     87.555  
13      -13.093835        -98978      0.106542    146.182  
14       13.093835         98978      0.106497   -146.239  
15        4.816002        269103      0.075031    149.249  
16       -4.816002       -269103      0.075030   -149.249  
17        2.100546        616982      0.114874    135.584  
18        0.337941       3834988      0.097371   -100.413  
19       -2.100547       -616982      0.114703   -135.554  
20       -0.338062      -3833617      0.097363    101.137  
21      -13.976711        -92726      0.052646     73.591  
22       13.976711         92726      0.051633    -73.535  
23        4.587860        282485      0.062891    157.242  
24       -4.587860       -282485      0.062891   -157.242  
25        2.361334        548842      0.073103    134.143  
26       -2.361334       -548842      0.073106   -134.132  
27        5.316474        243771      0.071482    168.324  
28        5.596725        231564      0.077239    104.389  
29       -5.316474       -243771      0.071482   -168.324  
30       -5.596725       -231564      0.077239   -104.389  
31        6.178106        209773      0.033163   -155.327  
32       -6.178106       -209773      0.033163    155.328  
33        1.891647        685117      0.058752     71.488  
34       -1.891650       -685116      0.059850    -70.705  
35        8.163185        158762      0.041085    128.291  
36       -8.163185       -158762      0.041085   -128.292  
37        3.815482        339669      0.066386    165.447  
38        8.567904        151262      0.066640     42.497  
39       -3.815482       -339669      0.066387   -165.445  
40       -8.567905       -151262      0.066639    -42.496  
41       -1.430111       -906223      0.085694    -49.790  
42       -1.183222      -1095315      0.080038     15.221  
43       -1.635218       -792555      0.073156   -155.506  
44      -13.771689        -94106      0.092299    -25.266
     */
/*f1=tf[2];f2=tf[7];f3=tf[14];f4=tf[10];w1=abs(za[2]);w2=abs(za[7]);w3=abs(za[14]);w4=abs(za[10]);*/
/* these are in arcsec/a - need to convert to Ma */
/*%subnafres[t_sr93_ip*1e6];*/
/*

*/

/* Frequencies from data*/

f1 = 3.2175636007322;
f2 = 13.6613294228515;
f3 = 13.1043002807929;
f4 = 10.4202281006986;
f5 = 9.87952708416141;
w1 = 0.356835893916081;
w2 = 0.106629598881501;
w3 = 0.124004190416502;
w4 = 0.146516892824212;
w5 = 0.0828652785594736;

f1 = 3.2266562897475;
f2 = 13.6480372568966;
f3 = 13.0938345214816;
f4 = 10.4154976442281;
w1 = 0.349414408361938;
w2 = 0.124519596404924;
w3 = 0.10649692423927;
w4 = 0.131428636797921;

/* Frequencies from ETP93 */
/*
f1 = 3.20756328363849;
f2 = 13.6197939752291;
f3 = 13.0905385672221;
f4 = 10.4137603288143;
f5 = 9.88051500103392;
w1 = 0.396315698570965;
w2 = 0.300111838274881;
w3 = 0.165459559203283;
w4 = 0.252273666616271;
w5 = 0.144207823155108;
*/
f1=%asec_to_Ma[f1];f2=%asec_to_Ma[f2];f3=%asec_to_Ma[f3];f4=%asec_to_Ma[f4];f5=%asec_to_Ma[f5];

wsum=w1+w2+w3+w4+w5;
/*filtwidth2=0.9;*/

w1dem=w1/wsum*%demod[f1_93,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];
w1dem=2*%avg[w1dem]-w1dem; /* flip 406 ky component */
w23dem=(w2+w3)/wsum*%demod[f23_93,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];
w45dem=(w4+w5)/wsum*%demod[f45_93,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];

w1dem=w1/wsum*%demod[f1,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];
w1dem=2*%avg[w1dem]-w1dem; /* flip 406 ky component */
w2dem=w2/wsum*%demod[f2,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];
w3dem=w3/wsum*%demod[f3,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];
w4dem=w4/wsum*%demod[f4,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];
w5dem=w5/wsum*%demod[f5,t_sr93_ip,sr93_pdem,-filtwidth2,filtwidth2];

sr93_prececc_dem=w1dem+w23dem+w45dem;
sr93_prececc_dem=w1dem+w2dem+w3dem+w4dem+w5dem;
sr93_prececc_dem=w1dem+w23dem+w45dem;
%filtreR[sr93_prececc_dem,t_sr93_ip,1,-final_filter_freq,final_filter_freq];
sr93_prececc_dem=_fi_x;


/*plot(t_sr93_ip,w1dem,opt);
replot(t_sr93_ip,w23dem,opt);
replot(t_sr93_ip,w45dem,opt);

plot(t_sr93_ip,w1dem,opt);
replot(t_sr93_ip,w2dem,opt);replot(t_sr93_ip,w3dem,opt);
replot(t_sr93_ip,w4dem,opt);replot(t_sr93_ip,w5dem,opt);*/
/*plot(t_sr93_ip,sr93_prececc_dem,opt);*/
/************END DEMDULATION FOR SR1993 *******************/

tmin=18;
tmax=29.0;

t_sr93_clip=select((t_sr93_ip <=tmax)&&(t_sr93_ip >=tmin),t_sr93_ip);
t_sr03_clip=select((t_sr03_ip <=tmax)&&(t_sr03_ip >=tmin),t_sr03_ip);
sr93_prececc_dem_clip=select((t_sr93_ip <= tmax)&&(t_sr93_ip >= tmin),sr93_prececc_dem);
sr03_prececc_dem_clip=select((t_sr03_ip <= tmax)&&(t_sr03_ip >= tmin),sr03_prececc_dem);


plot(t03_ip,etp03_prececc_dem,"title 'etp03' w l");
replot(t93_ip,etp93_prececc_dem,"title 'etp93' w l");
replot(t_sr03_clip,sr03_prececc_dem_clip+1,"title 'sr03' w l");
replot(t_sr93_clip,sr93_prececc_dem_clip+1,"title 'sr93' w l");

write("geol_prececc_93_03_res.txt",t_sr03_clip,sr93_prececc_dem_clip,sr03_prececc_dem_clip);

write("astr_prec_93_03_prececcdem.txt",t03_ip,etp93_pdem,etp03_pdem);
write("geol_prec_93_03_prececcdem.txt",t_sr03_ip,sr93_pdem,sr03_pdem);

};









macro eccentricity_demodulation{
/* demodulate ETP target and Leg 199 data at frequencies
identified with subnafr */
vnumR etp03_edem,etp93_edem,la199_edem;
vnumR etp03ecc,etp93ecc,la199ecc;
vnumR la199ecc_hil,etp93ecc_hil,etp03ecc_hil;
vnumR la199_edem_frac,etp93_edem_frac,etp03_edem_frac;
vnumR filtwidth;

filtwidth=0.6;

etp03_edem=((%demod[8.068420203,t03_ip,etp03_ip,-filtwidth,filtwidth]+%demod[10.53596451,t03_ip,etp03_ip,-filtwidth,filtwidth]))/2;
etp93_edem=((%demod[8.033741715,t93_ip,etp93_ip,-filtwidth,filtwidth]+%demod[10.10325527,t93_ip,etp93_ip,-filtwidth,filtwidth]))/2;
la199_edem=((%demod[8.042593576,t_199_ip,dat_199_ip_sm,-filtwidth,filtwidth]+%demod[10.53529852,t_199_ip,dat_199_ip_sm,-filtwidth,filtwidth]))/2;

%filtreR[etp03_ip,t03_ip,1,6,12];
etp03ecc=_fi_x;
%filtreR[etp93_ip,t93_ip,1,6,12];
etp93ecc=_fi_x;
%filtreR[dat_199_ip_sm,t_199_ip,1,6,12];
la199ecc=_fi_x;

la199ecc_hil=abs(%hilbert[la199ecc]);
la199ecc_hil=%smooth2[1500,la199ecc_hil];
la199_edem_frac=la199_edem/la199ecc_hil;

etp93ecc_hil=abs(%hilbert[etp93ecc]);
etp93ecc_hil=%smooth2[1500,etp93ecc_hil];
etp93_edem_frac=etp93_edem/etp93ecc_hil;

etp03ecc_hil=abs(%hilbert[etp03ecc]);
etp03ecc_hil=%smooth2[1500,etp03ecc_hil];
etp03_edem_frac=etp03_edem/etp03ecc_hil;

plot(t_199_ip,la199_edem_frac,opt);
replot(t03_ip,etp03_edem_frac,opt);
replot(t93_ip,etp93_edem_frac,opt);
write("astr_ecc_93_03_res.txt",t03_ip,etp93_edem_frac,etp03_edem_frac);
write("geol_ecc_03_res.txt",t_199_ip,la199_edem_frac);

write("geol_ecc_03_edem.txt",t_199_ip,la199_edem);


};



/*********************** MACROS ************************/
macro subnafr[_tt,_xx]{
/* analyse en frequence */
vnumR tf,trx,try;
vnumC za;

_naf_dtour=360*3600;

_yy=_xx*0;

_step = (_tt[size(_tt)]-_tt[1])/(size(_xx)-1);

naftab(_xx,_yy,za,tf,size(_xx),_step,_tt[1],NTERM,trx,try);
P=360*3600/tf;
ang=arg(za)/PI*180;
writes("%15.6f  %12.0f  %12.6f %10.3f  \n",tf,P,abs(za),ang); 
};

macro subnaf2[_tt,_xx,_yy]{
/* analyse en frequence */
vnumR tf,trx,try;
vnumC za;

_naf_dtour=360*3600;

_step = (_tt[size(_tt)]-_tt[1])/(size(_xx)-1);

naftab(_xx,_yy,za,tf,size(_xx),_step,_tt[1],100,trx,try);
P=360*3600/tf;
ang=arg(za)/PI*180;
writes("%15.6f  %12.0f  %12.6f %10.3f  \n",tf,P,abs(za),ang); 
};

macro subnafres[_tt]{
/* analyse en frequence */

naftab(trx,try,za,tf,size(trx),_step,_tt[1],100,trx,try);
P=360*3600/tf;
ang=arg(za)/PI*180;
writes("%15.6f  %12.0f  %12.6f %10.3f  \n",tf,P,abs(za),ang); 
};

macro plotfftR[_xx,_tt,_dtour,_fmin,_fmax,_IW] {
/* trace de la fft d'un signal reel */

vnumR _ftf;
vnumC _fza;

_yy=_xx*0;

_step = (_tt[2]-_tt[1]);
fft(_xx,_yy,_fza,_ftf,_step,_tt[1],_dtour,_IW);

strange=str(_fmin)+":"+str(_fmax);
gnuplot;
set xrange [@strange@]
end;
plot(_ftf,abs(_fza)," w l");
};


macro filtreR[_xx,_tt,_dtour,_fmin,_fmax] {
/* filtre un signal reel _xx    jxl 7/2/2002*/
vnumR _ftf;
vnumC _fza,_fi_fza;
vnumR _ix,_iy;
vnumR _fi_x,_fi_y;

_step = _tt[2]-_tt[1];
_IW=0;

_yy=_xx*0;
/* analyse de Fourier */
fft(_xx,_yy,_fza,_ftf,_step,_tt[1],_dtour,_IW);

/* filtre */
_fi_fza=?((abs(_ftf )> _fmin) && (abs(_ftf) < _fmax)):_fza:0;

/* analyse de Fourier inverse*/
ifft(_fi_fza,_fi_x,_fi_y,_step,_tt[1]);

};

macro filtreC[_cc,_tt,_dtour,_fmin,_fmax] {
/* filtre un signal complexe _cc    hp 6/2/2003*/
vnumR _ftf;
vnumC _fza,_fi_fza;
vnumR _ix,_iy;
vnumR _fi_x,_fi_y;
vnumR _xx,_yy;

_step = _tt[2]-_tt[1];
_IW=0;
_xx=real(_cc);
_yy=imag(_cc);
/* analyse de Fourier */
fft(_xx,_yy,_fza,_ftf,_step,_tt[1],_dtour,_IW);

/* filtre */
_fi_fza=?((abs(_ftf )> _fmin) && (abs(_ftf) < _fmax)):_fza:0;

/* analyse de Fourier inverse*/
ifft(_fi_fza,_fi_x,_fi_y,_step,_tt[1]);

};

macro smooth2[_pas,_xx]{
/*********************************************************/
/* effectue un lissage simple par moyenne sur -_pas,_pas */
/*********************************************************/
_sxx=_xx;
_siz=size(_xx);
for n=1 to _pas{
   _sxx[n ] = sum(_xx[1:n+_pas])/(n+_pas) $};
   
for n = _pas+1 to _siz -_pas {
   _sxx[n ] = sum(_xx[n-_pas:n+_pas])/(2*_pas+1) $};

for n = _siz -_pas+1 to _siz {
   _sxx[n ] = sum(_xx[n-_pas:_siz])/(_siz-n+_pas+1) $};

return(_sxx);
};

macro smooth3[_pas,_xx]{
/*********************************************************/
/* effectue un lissage simple par moyenne sur -_pas,_pas 
   nouvelle version pour faire mieux les bords           */
/*********************************************************/
_siz=size(_xx);
/* addition de donnees des 2 cotes */
_AT =1,_pas;
_DT0=_xx[1:_pas];
%least[_AT,_DT0];
_ATn = -_pas+1,0;

_AT0=a*_ATn + b;

_DT1=_xx[_siz-_pas+1:_siz];
%least[_ATn,_DT1];
_AT1=a*_AT + b;

_gxx=vnumR[_AT0:_xx:_AT1];

_sgxx=%smooth2[_pas,_gxx];
_sxx=_sgxx[_pas+1:_siz-_pas];
return(_sxx);
};

macro hilbert[_xx]{
/****************************/
/* HILBERT TRANSFORM */
/* return complex hilbert transform */
/* assumes equally spaced xx */

_yy=_xx*0;


vnumC z;
vnumR freq;
vnumR n;
vnumC h,zrh,hilbertres;

n=2^int(log(size(_xx))/log(2));


/* step 1: compute fft of _xx */
fft(_xx,_yy,z,freq,1,0,2*pi,0);

/* step 2: compute h such that
It creates a vector hwhose elements h(i) have the values: 
1 for iÊ=Ê DC part
2 for iÊ=Ê pos freq
0 for iÊ=Ê neg freq
*/
resize(h,n+1,0); 
/* index for DC part*/
ind_f0=imin(abs(freq));
h[ind_f0]=1;
/* set h=2 for pos frequencies */
h[(ind_f0+1):(n+1)]=2;
/* all others are zero already */
/*It calculates the element-wise product of xand h.*/
zrh=z*h;
/*It calculates the inverse FFT of xrh and returns the first nelements of the result */
ifft(zrh,xout,yout,1,0);
hilbertres=_xx+i*yout;
return(hilbertres);
};

macro derivee[T]{
T1=T[2:size(T)];
T2=T[1:size(T)-1];
DT=T1-T2;
return(DT);
};


macro stdev[x]{
vnumR xavg,xsqsum,xsumsq,N,xstdev;
N=size(x);
xavg=sum(x)/N;
xsqsum=sum(x^2);
xsumsq=sum(x)^2;
xstdev=sqrt((N*xsqsum-xsumsq)/(N*(N-1)));
return(xstdev);
};

macro avg[x]{
return(sum(x)/size(x));
};

macro asec_to_Ma[x]{
return(1/(360*3600/x/1e6));
};


macro unwrap_phase[_ph]{
vnumR unwrapped_phase;
phPrev=0;
phCum=0;
phNew=0;
delta=0;
length=size(_ph); 
resize(unwrapped_phase,length);
for j=1 to length
{ 
   phNew = _ph[j]$ 
   delta = phNew-phPrev$
   phPrev = phNew$
   if(delta > pi) then 
     {delta = delta - 2*pi$}
   else 
     {if(delta < -pi) then {delta = delta + 2*pi$};
   };
   phCur = delta + phCum$ 
   phCum = phCur$
   unwrapped_phase[j] = phCur$
 };
 return(unwrapped_phase);
}; 



/*********************** END OF MACROS ************************/
