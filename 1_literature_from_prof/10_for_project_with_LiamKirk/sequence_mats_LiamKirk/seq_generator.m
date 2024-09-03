% dim_size  : size of the lattice, e.g. for 8^4 then pass a value of 8 here
% config_nr : for a particular lattice size, the id of the configuration

% example use : seq_generator(4,0,0.01,4)

function [] = seq_generator(dim_size, config_nr, delta_seq, seq_length)

  % only first config at the moment
  if config_nr~=0
    error("Only the first config is used at the moment");
  end

  in_dir = "../QCDconfigs/";
  out_dir = "../";

  switch dim_size
   case 4
      beta = 3.55;
      kappa = 0.137;
      %csw = 1.8248654;
      csw = 0.0;
      chem_pot = 0.3;
      if config_nr<0 || config_nr>19
          error("For 4^4, only configs going from 0 to 19");
      end
   case 8
      beta = 3.55;
      kappa = 0.137;
      %csw = 1.8248654;
      csw = 0.0;
      chem_pot = 0.3;
      if config_nr<0 || config_nr>19
          error("For 8^4, only configs going from 0 to 19");
      end
   case 16
      beta = 3.55;
      kappa = 0.137;
      %csw = 1.8248654;
      csw = 0.0;
      chem_pot = 0.3;
      if config_nr<0 || config_nr>19
          error("For 16^4, only configs going from 0 to 19");
      end
   case 32
      beta = 3.55;
      kappa = 0.137;
      %csw = 1.8248654;
      csw = 0.0;
      chem_pot = 0.3;
      if config_nr<0 || config_nr>19
          error("For 32^4, only configs going from 0 to 19");
      end
  otherwise
      error("The value of dim size can be only 16 or 32 for now");
  end

  out_dir = strcat(out_dir,string(dim_size));
  out_dir = strcat(out_dir,"to4/");

  m0 = 1/(2*kappa)-4;
  
  file_name = "periodic_L";
  file_name = strcat(file_name,string(dim_size));
  file_name = strcat(file_name,"_b");
  file_name = strcat(file_name,string(beta));
  file_name = strcat(file_name,"_k");
  file_name = strcat(file_name,string(kappa));
  file_name = strcat(file_name,"n");
  file_name = strcat(file_name,string(config_nr));

  fprintf("Converting : %s, located at %s\n",file_name,in_dir);

  lat = 0;
  %delta_seq = 0.01;
  for ix=1:seq_length
    [D,lat] = wilson_dirac_from_cnfg(strcat(in_dir,file_name),m0,csw,chem_pot,lat,ix,delta_seq);

    out_name = strcat(file_name,"_");
    out_name = strcat(out_name,string(ix));
    out_name = strcat(out_name,".mat");

    fprintf("Saving to : %s\n",strcat(out_dir,out_name));
    save(strcat(out_dir,out_name),'D','-v7.3');
  end

end
