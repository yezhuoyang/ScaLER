from LERcalc.LERcalculator import *










distance=3
circuit=CliffordCircuit(2)
circuit.set_error_rate(0.00001)
stim_circuit=stim.Circuit.generated("surface_code:rotated_memory_z",rounds=distance*3,distance=distance).flattened()
stim_circuit=rewrite_stim_code(str(stim_circuit))
circuit.set_stim_str(stim_circuit)
circuit.compile_from_stim_circuit_str(stim_circuit)           
new_stim_circuit=circuit.get_stim_circuit()      


# Configure a decoder using the circuit.
detector_error_model = new_stim_circuit.detector_error_model(decompose_errors=False)
matcher = pymatching.Matching.from_detector_error_model(detector_error_model)



result=return_samples(str(new_stim_circuit),4,1000000)
