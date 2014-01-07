function data = readClojureSimFile(filename)
% data = readClojureSimFile(filename)
% Reads the binary formatted file output by the clojure retina simulation code. Yields a struct with
% members labels, mosaic, and R representing the cone types, cone positions, and correlation matrix,
% respectively.
   f = fopen(filename, 'r');
   try
       n = fread(f, 1, 'int32', 0, 'ieee-be');
       % the next n integers are cone labels...
       data.labels = fread(f, n, 'int32', 0, 'ieee-be');
       % then, there are 2*n 4-byte floats for cone locations
       data.mosaic = fread(f, [2 n], 'float32', 0, 'ieee-be');
       % then, there's the correlation matrix
       data.R = fread(f, [n n], 'float32', 0, 'ieee-be');
       % last, there's the embedding (if present)
       data.em = fread(f, [n 3], 'float32', 0, 'ieee-be');
   catch e
       fclose(f);
       rethrow(e);
   end
   fclose(f);
end
