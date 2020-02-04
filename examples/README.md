Explanation of Example Files
============================
David Brown

9 Dec. 2019

Each step has a "sample_input" containing all needed inputs and a
"sample_output" with the expected output.  For example, `chum` produces
`test_phasemaker_output.dat` in `CHUM/sample_output/` and this file is identical to `CRAB/sample_input/test_phasemaker_output.dat`.

| Step | Code       | Input files | Output files |
| ---- | ---------- | ----------- | -------------|
|  1   | chum       | test_phasemaker_input.dat | test_phasemaker_output_oscar.dat |
|  2   | crab       | test_crab_input.dat  test_phasemaker_output_oscar.dat|  test_crab_output_correlation.dat  test_crab_output_source.dat|   
|  2'  | converts2c |             |              |
|  3   | scplot     |             |              |
|  4   | bfplot     |             |              |
|  5   | diver      |             |              |
|  6   | scplot     |             |              |
|  7   | converts2c |             |              |
