%============================================================================
%   swmm.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01          (Build 2016/01/08 by L. Rossma)
%   Author:     L. Rossma
%   Modified by:Vinh Ngoc Tran
%
%   Simulation of swmm5
%============================================================================
classdef swmm < handle
  properties (Constant)
    JUNCTION = 0;
    SUBCATCH = 1;
    NODE = 2;
    LINK = 3;
    STORAGE = 4;
    ORIFICE = 414;
    US = 0;
    SI = 1;
    DIMENTIONLESS = 0;
    DEPTH = 200;
    VOLUME = 201;
    FLOW = 202;
    SETTING = 203;
    FROUDE = 204;
    INFLOW = 205;
    FLOODING = 206;
    PRECIPITATION = 207;
    RUNOFF = 208;
    MAX_AREA = 209;
    CAPACITY = -200;
    NO_REPORT = 0;
    WRITE_REPORT = 1;
    INVERT = 400;
    DEPTH_SIZE = 401;
    STORAGE_A = 402;
    STORAGE_B = 403;
    STORAGE_C = 404;
    LENGTH = 405;
    ROUGHNESS = 406;
    IN_OFFSET = 407;
    OUT_OFFSET = 408;
    AREA = 409;
    IMPERV = 410;
    WIDTH = 411;
    SLOPE = 412;
    OUTLET = 413;
    FROM_NODE = 415;
    TO_NODE = 416;
    NONE = -400;
  end
  properties
    elapsed_time;
    timePtr;
    is_initialized = false;
  end
  properties (Hidden = true)
    % Error codes
    ERROR_PATH = -300;
    ERROR_ATR = -299;
    ERROR_TYPE = -298;
    ERROR_NFOUND = -297;
    ERROR_INCOHERENT = -296;
    ERROR_IS_NUMERIC = -295;
    % Error messages - Exceptions
    ERROR_MSG_NFOUND = MException('AttributeError:Check_ID', ...
            'Error: Object not found');
    ERROR_MSG_TYPE = MException('AttributeError:Check_TYPE', ...
            'Error: Type of object not compatible');
    ERROR_MSG_ATR = MException('AttributeError:Check_ATRBT', ...
            'Error: Attribute not compatible');
    ERROR_MSG_PATH = MException('AttributeError:Check_FILE_PATH', ...
            'Error: Incorrect file path');
    ERROR_MSG_INCOHERENT = MException('TypeError:Check_PARAMETERS', ...
            'Error: Incoherent parameter');
    ERROR_MSG_SYSTEM = MException('systemError:SystemFailure', ...
            'Error: The system failed - files must be closed');
    ERROR_MSG_IS_NUMERIC = MException('AttributeError:NotNumeric', ...
            'Error: This function just handle numerical attributes');
  end
  %%
  methods
  %%
  function open(obj, input_file)
    if ~(libisloaded('swmm5'))
      loadlibrary('swmm5');
    end

    rpt_file = strrep(lower(input_file), '.mdf', '.rpt');
    out_file = strrep(lower(input_file), '.mdf', '.out');
    error = calllib('swmm5','swmm_open',input_file, rpt_file, out_file);
    if error ~= 0
      if (libisloaded('swmm5'))
        unloadlibrary swmm5;
      end
      throw(obj.ERROR_MSG_PATH);
    end
  end
  %%
  function start(obj, write_report)
    if ~ismember(write_report, [obj.NO_REPORT, obj.WRITE_REPORT])
      if libisloaded('swmm5')
        unloadlibrary swmm5;
      end
      throw(obj.ERROR_MSG_INCOHERENT);
    end
    if ~(libisloaded('swmm5'))
      loadlibrary('swmm5');
    end
    error = calllib('swmm5','swmm_start', write_report);
    obj.elapsed_time = 1e-6;
    obj.timePtr = libpointer('doublePtr', obj.elapsed_time);

    if error ~= 0
      if libisloaded('swmm5')
        unloadlibrary swmm5;
      end
      exception = MException('SystemFailure:CheckErrorCode',...
      sprintf('Error %d ocurred', error));
      throw(exception);
    end

    obj.is_initialized = true;
  end
  %%
  function time = run_step(obj)
    if ~(libisloaded('swmm5'))
      loadlibrary('swmm5');
    end
    error = calllib('swmm5','swmm_step', obj.timePtr);
    time  = obj.timePtr.value*24;

    if error ~= 0
      exception = MException('SystemFailure:CheckErrorCode',...
      sprintf('Error %d ocurred at time %.2f hours', error, time));
      if libisloaded('swmm5')
        unloadlibrary swmm5;
      end
      throw(exception);
    end
  end
  %%
  function duration = end_sim(obj)
    if ~(libisloaded('swmm5'))
      loadlibrary('swmm5');
    end

    error = calllib('swmm5','swmm_end');

    if error ~= 0
      exception = MException('SystemFailure:CheckErrorCode',...
      sprintf('Error %d: The simulation can not be ended', error));
      if libisloaded('swmm5')
        unloadlibrary swmm5;
      end
      throw(exception);
    end
    duration = 0;
    obj.is_initialized = false;
  end
  %%
  function report(obj)
    if ~(libisloaded('swmm5'))
      loadlibrary('swmm5');
    end

    error = calllib('swmm5','swmm_report');

    if error ~= 0
      exception = MException('SystemFailure:CheckErrorCode',...
      sprintf('Error %d: The report file could not be written correctly', error));
      if libisloaded('swmm5')
        unloadlibrary swmm5;
      end
      throw(exception);
    end
  end
  %%
  function close(obj)
    if ~(libisloaded('swmm5'))
      loadlibrary('swmm5');
    end

    error = calllib('swmm5','swmm_close');

    if error ~= 0
      exception = MException('SystemFailure:CheckErrorCode',...
      sprintf('Error %d: The file can not be closed correctly', error));
      if libisloaded('swmm5')
        unloadlibrary swmm5;
      end
      throw(exception);
    end
    if libisloaded('swmm5')
      unloadlibrary swmm5;
    end
  end
  %%
  function errors = get_mass_bal_error(obj)
    runOffErr = single(0);
    flowErr = single(0);
    qualErr = single(0);
    ptrRunoff = libpointer('singlePtr', runOffErr);
    ptrFlow = libpointer('singlePtr', flowErr);
    ptrQual = libpointer('singlePtr', qualErr);

    if ~(libisloaded('swmm5'))
      loadlibrary('swmm5');
    end

    error = calllib('swmm5','swmm_getMassBalErr', ptrRunoff, ptrFlow, ptrQual);
    if error ~= 0
      exception = MException('SystemFailure:CheckErrorCode',...
      sprintf('Error %d: The errors can not be retrieved', error));
      if libisloaded('swmm5')
        unloadlibrary swmm5;
      end
      throw(exception);
    end

    runoff = ptrRunoff.value;
    flow = ptrFlow.value;
    qual = ptrQual.value;
    errors = [runoff, flow, qual];
  end
  %%
  function bool_ans = is_over(obj)
    bool_ans = obj.timePtr.value == 0;
  end
  %%
  function values = step_get(obj, list_ids, attribute, unit_system)
    if isa(list_ids, 'char')
      values = obj.get(list_ids, attribute, unit_system);
    else
      values = zeros(1, length(list_ids));
      for i=1 : length(list_ids)
        values(i) = obj.get(list_ids{i}, attribute, unit_system);
      end
    end
  end
  %%
  function value = get_from_input(obj, input_file, object_id, attribute)
 
  if ~ismember(attribute, [obj.INVERT, obj.DEPTH_SIZE, obj.STORAGE_A, ...
  obj.STORAGE_B, obj.STORAGE_C, obj.LENGTH, obj.ROUGHNESS, ...
  obj.IN_OFFSET, obj.OUT_OFFSET, obj.AREA, obj.IMPERV,  obj.WIDTH, obj.SLOPE])
  throw(obj.ERROR_MSG_INCOHERENT);
  end

  if ~(libisloaded('swmm5'))
  loadlibrary('swmm5');
  end

  value = calllib('swmm5','swmm_get_from_input', input_file, object_id, attribute);

  if value == obj.ERROR_NFOUND
  throw(obj.ERROR_MSG_NFOUND);
  elseif value == obj.ERROR_TYPE
  throw(obj.ERROR_MSG_TYPE);
  elseif value == obj.ERROR_ATR
  throw(obj.ERROR_MSG_ATR);
  elseif value == obj.ERROR_PATH
  throw(obj.ERROR_MSG_PATH);
  end
  end
  %%
  function [id_list, value_list] = get_all(obj, input_file, object_type, attribute)
  
  if (attribute == obj.NONE)
  attribute = -1;
  end

  if (attribute == obj.MAX_AREA)
  if object_type ~= obj.LINK
  throw(obj.ERROR_MSG_ATR);
  end
  id_list = obj.get_all(input_file, object_type, obj.NONE);
  obj.initialize(input_file);
  for i = 1 : length(id_list)
  value_list(i) = obj.get(id_list{i}, obj.MAX_AREA, obj.SI);
  end
  obj.finish;
  return;
  end

  if attribute ~= -1
  if ~ismember(attribute, [obj.INVERT, obj.DEPTH_SIZE, obj.STORAGE_A, ...
  obj.STORAGE_B, obj.STORAGE_C, obj.LENGTH, obj.ROUGHNESS, ...
  obj.IN_OFFSET, obj.OUT_OFFSET, obj.AREA, obj.IMPERV,  obj.WIDTH, obj.SLOPE, obj.OUTLET, ...
  obj.FROM_NODE, obj.TO_NODE])
  throw(obj.ERROR_MSG_INCOHERENT);
  end
  end
  if ~ismember(object_type, [obj.JUNCTION, obj.SUBCATCH, obj.LINK, obj.STORAGE, obj.ORIFICE, obj.NODE])
  throw(obj.ERROR_MSG_INCOHERENT);
  end

  if ~(libisloaded('swmm5'))
  loadlibrary('swmm5');
  end

  error = calllib('swmm5','swmm_save_all', input_file, object_type, attribute);

  if(error == obj.ERROR_PATH)
  throw(obj.ERROR_MSG_PATH);
  elseif (error == obj.ERROR_NFOUND)
  delete('info.dat');
  throw(obj.ERROR_MSG_NFOUND);
  elseif (error == obj.ERROR_ATR)
  delete('info.dat');
  throw(obj.ERROR_MSG_ATR);
  elseif (error == obj.ERROR_TYPE)
  delete('info.dat');
  throw(obj.ERROR_MSG_TYPE);
  end

  if(attribute == -1)
  IS_VECTOR = true;
  else
  IS_VECTOR = false;
  end

  file = fopen('info.dat','r');
  i = 1;
  if IS_VECTOR
  while ~feof(file)
  id_list{i} = fgetl(file);
  i = i+1;
  end
  if attribute ~= obj.OUTLET
  value_list = [];
  else
  value_list = {};
  end
  else
  while ~feof(file)
  line = fgetl(file);
  if line ~= -1
      if ~ismember(attribute, [obj.OUTLET obj.FROM_NODE obj.TO_NODE])
          key_value     = textscan(line, '%s %f');
          id_list{i}    = char(key_value{1});
          value_list(i,:) = key_value{2};
      else
          key_value     = textscan(line, '%s %s');
          id_list{i}    = char(key_value{1});
          value_list{i} = char(key_value{2});
      end

      i = i+1;
  end
  end

  end
  fclose(file);
  delete('info.dat');
  end
  %%
  function modify_input(obj, input_file, object_id, attribute, value)
  
  if ~(libisloaded('swmm5'))
  loadlibrary('swmm5');
  end

  error = calllib('swmm5','swmm_modify_input', input_file, object_id, attribute, value);

  if error == obj.ERROR_NFOUND
  throw(obj.ERROR_MSG_NFOUND);
  elseif error == obj.ERROR_TYPE
  throw(obj.ERROR_MSG_TYPE);
  elseif error == obj.ERROR_ATR
  throw(obj.ERROR_MSG_ATR);
  elseif error == obj.ERROR_PATH
  throw(obj.ERROR_MSG_PATH);
  end
  end
  %%
  function modify_settings(obj, orifices_ids, new_settings)
  
  for i=1:length(orifices_ids)
  obj.modify_setting(orifices_ids{i}, new_settings(i));
  end
  end
  %%

  function [errors, duration] = run_simulation(obj, input_file)
  
  obj.initialize(input_file);
  while ~obj.is_over
  obj.run_step;
  end
  [errors, duration] = obj.finish;
  end
  function initialize(obj, input_file)
  obj.open(input_file);
  obj.start(obj.WRITE_REPORT);
  end
  %%
  function save_results(obj)
  error = calllib('swmm5','swmm_save_results');
  end
  %%
  function [errors, duration] = finish(obj)
  duration = obj.end_sim;
  errors = obj.get_mass_bal_error;
  obj.report;
  if exist('Nodes') == 7
  rmdir('Nodes', 's');
  end
  if exist('Links') == 7
  rmdir('Links', 's');
  end
  if exist('Subcatchments') == 7
  rmdir('Subcatchments', 's');
  end
  if exist('Time') == 7
  rmdir('Time', 's');
  end
  obj.save_results;
  obj.close;
  end
  %%
  function [time, result] = read_results(obj, object_id, object_type, attribute)
  if object_type == obj.LINK
  folder = 'Links/';
  compatible = [obj.FLOW, obj.DEPTH, obj.VOLUME, obj.CAPACITY];
  columns = [1, 3, 4, 5];
  elseif object_type == obj.NODE
  folder = 'Nodes/';
  compatible = [obj.INFLOW, obj.FLOODING, obj.DEPTH, obj.VOLUME];
  columns = [1, 2, 3, 5];
  elseif object_type == obj.SUBCATCH
  folder = 'Subcatchments/';
  compatible = [obj.PRECIPITATION, obj.RUNOFF];
  columns = [1, 3];
  else
  throw(obj.ERROR_MSG_TYPE);
  end
  [a, position] = ismember(attribute, compatible);
  if a ~= 1
  throw(obj.ERROR_MSG_ATR);
  end

  if ~iscell(object_id)
  object_id = {object_id};
  end
  result = [];
  for i=1 : length(object_id)
  path = strcat(folder, object_id{i}, '.csv');
  if exist(path) ~= 2
  throw(obj.ERROR_MSG_NFOUND);
  end
  data = csvread(path);
  result(:,i) = [0;data(:,columns(position))];
  end
  info = csvread('Time/time.txt');
  time = zeros(info(2),1);
  for i=2 : info(2)+1
  time(i) = info(1)/3600 + time(i-1);
  end
  end
  %%
  function modify_setting(obj, orifice_id, new_setting)
  if ~(libisloaded('swmm5'))
  loadlibrary('swmm5');
  end

  error = calllib('swmm5','swmm_modify_setting', orifice_id, new_setting, 0);
  if error == obj.ERROR_INCOHERENT
  throw(obj.ERROR_MSG_INCOHERENT);
  elseif error == obj.ERROR_NFOUND
  throw(obj.ERROR_MSG_NFOUND);
  end
  end
  %%
  function value = get(obj, object_id, attribute, unit_system)
  if ~ismember(attribute, [obj.DEPTH, obj.VOLUME, obj.FLOW, obj.SETTING, obj.FROUDE, obj.INFLOW, obj.FLOODING, ...
      obj.PRECIPITATION, obj.RUNOFF, obj.MAX_AREA])
  throw(obj.ERROR_MSG_INCOHERENT);
  elseif ~ismember(unit_system, [obj.SI, obj.US, obj.DIMENTIONLESS])
  throw(obj.ERROR_MSG_INCOHERENT);
  end
  if ~(libisloaded('swmm5'))
  loadlibrary('swmm5');
  end
  value = calllib('swmm5','swmm_get', object_id, attribute, unit_system);
  if (value == obj.ERROR_NFOUND)
  throw(obj.ERROR_MSG_NFOUND);
  elseif (value == obj.ERROR_TYPE)
  throw(obj.ERROR_MSG_TYPE);
  elseif (value == obj.ERROR_ATR)
  throw(obj.ERROR_MSG_ATR);
  end
  end
  end
end
