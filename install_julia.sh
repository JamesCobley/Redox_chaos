# Install Julia in Colab
curl -fsSL https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.3-linux-x86_64.tar.gz -o julia.tar.gz
tar -xzf julia.tar.gz
mv julia-1.9.3 /usr/local/julia
ln -s /usr/local/julia/bin/julia /usr/local/bin/julia
