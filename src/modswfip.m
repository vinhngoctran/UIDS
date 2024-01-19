%============================================================================
%   modswfip.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program reads inputs for Soil water fluxes model
%============================================================================

function modsfl_inf = modswfip(filename)

% This file is for only read lumped information with only one soil layer for soil water flux model.
% it should be modified if apply for distributed area

%% Read textfile
opts = delimitedTextImportOptions("NumVariables", 20);
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["nlayer", "VarName2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20"];
opts.SelectedVariableNames = ["nlayer", "VarName2"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["nlayer", "VarName2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["nlayer", "VarName2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20"], "EmptyFieldRule", "auto");
SOILFLUX1 = readmatrix(filename, opts);
clear opts;

%% Extract input
modsfl_inf.nlayer = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["nlayer"]),2));            
modsfl_inf.thick = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["thick"]),2));                         
modsfl_inf.stonef = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["stonef"]),2));          
modsfl_inf.psicr = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["psicr"]),2));           
modsfl_inf.par(1,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par1"]),2));
modsfl_inf.par(2,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par2"]),2)); 
modsfl_inf.par(3,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par3"]),2)); 
modsfl_inf.par(4,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par4"]),2)); 
modsfl_inf.par(5,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par5"]),2)); 
modsfl_inf.par(6,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par6"]),2)); 
modsfl_inf.par(7,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par7"]),2)); 
modsfl_inf.par(8,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par8"]),2)); 
modsfl_inf.par(9,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par9"]),2)); 
modsfl_inf.par(10,1) = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["par10"]),2)); 
modsfl_inf.ilayer = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["ilayer"]),2));

modsfl_inf.infexp = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["infexp"]),2));           
modsfl_inf.qlayer = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["qlayer"]),2));            
modsfl_inf.dispc = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["dispc"]),2));           
modsfl_inf.height = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["height"]),2));                                     
modsfl_inf.lai = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["lai"]),2));                                    
modsfl_inf.sai = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["sai"]),2));            
modsfl_inf.snow = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["snow"]),2));              
modsfl_inf.snoden = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["snoden"]),2));           
modsfl_inf.mxrtln = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["mxrtln"]),2));           
modsfl_inf.mxkpl = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["mxkpl"]),2));  

modsfl_inf.densef = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["densef"]),2));        
modsfl_inf.frelden = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["frelden"]),2));            
modsfl_inf.tini = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["tini"]),2));               
modsfl_inf.age = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["age"]),2));               
modsfl_inf.rgroper = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["rgroper"]),2));           
modsfl_inf.inirdep = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["inirdep"]),2));         
modsfl_inf.inirlen = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["inirlen"]),2));          
modsfl_inf.rtrad = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["rtrad"]),2));           
modsfl_inf.fxylem = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["fxylem"]),2));  

modsfl_inf.qffc = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["qffc"]),2));            
modsfl_inf.qfpar = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["qfpar"]),2));        
modsfl_inf.dslope = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["dslope"]),2));        
modsfl_inf.length = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["length"]),2));           
modsfl_inf.drain = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["drain"]),2));         
modsfl_inf.gsc = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["gsc"]),2));           
modsfl_inf.gsp = str2num(SOILFLUX1(find(SOILFLUX1(:,1) == ["gsp"]),2));  

end