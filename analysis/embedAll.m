function count = embedAll(where)
% count = embedAll(where)
% This matlab function will read all of the .bin files in the given directory and run the
% mdscale function on their correlation matrices.  The results are placed in a large cell
% array em such that for em{1,i} == filename of ith file and em{2,i} == embedding of ith
% file.  This is then written to an identical file in the current directory.
% By default, if where is not given, it is assumed to be '.'.  The count value returned is
% the number of embeddings successfully processed.

  if nargin == 0, where = '.'; end;
  % process filenames first
  tmp = dir(where);
  flnames = {};
  for flno = 1:numel(tmp)
      flnm = tmp(flno).name;
      if flnm(1) == '.' || ~isequal(flnm(end-3:end), '.bin'), continue; end;
      flnames{end+1} = flnm;
  end
  % read and process the files...
  count = 0;
  for i = 1:numel(flnames)
      fprintf('%s...\n', flnames{i});
      q = readClojureSimFile([where '/' flnames{i}]);
      %tmp = find(q.R < 1);
      %mn = min(q.R(tmp));
      %mx = max(q.R(tmp));
      %q.R(tmp) = 0.925 * (q.R(tmp) - mn) / (mx - mn) + 0.05;
      q.R(find(q.R > 1 & q.R < 1.01)) = 1;
      q.em = mdscale(-log((q.R + 1) / 2), 3);
      writeClojureSimFile(['./' flnames{i}], q);
      count = count + 1;
  end
end
