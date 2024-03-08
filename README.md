# SimplicialExtended
Library for simplex complex & hypergraph.


## Usage

```
cmake -S . -B build
cd build
cmake --build .
```

This generate unit test and python module (simpl). 

## Generating clique Complex
For testing you can generate letter_image dataset using

```
python3 ./testing/generate_image_dataset.py.
```

Then build project using cmake:
```
cmake -S . -B build
cd build
cmake --build .
```

Then launch ./simplex_test