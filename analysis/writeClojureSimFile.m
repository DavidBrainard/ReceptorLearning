function writeClojureSimFile(filename, dat)
% data = readClojureSimFile(filename, dat)
% Writes a binary formatted clojure simulation file from the clojure simulation data in dat.
   f = fopen(filename, 'w');
   try
       fwrite(f, numel(dat.labels), 'int32', 0, 'ieee-be');
       % the next n integers are cone labels...
       fwrite(f, dat.labels, 'int32', 0, 'ieee-be');
       % then, there are 2*n 4-byte floats for cone locations
       fwrite(f, dat.mosaic, 'float32', 0, 'ieee-be');
       % then, there's the correlation matrix
       fwrite(f, dat.R, 'float32', 0, 'ieee-be');
       % finally, the embedding
       if ~isempty(dat.em)
           fwrite(f, dat.em, 'float32', 0, 'ieee-be');
       end
   catch e
       fclose(f);
       rethrow(e);
   end
   fclose(f);
end
