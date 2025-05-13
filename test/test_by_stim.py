import stim
import numpy as np


from QEPG.QEPG import return_samples,return_samples_many_weights,return_detector_matrix




oneQGate_ = ["H", "P", "X", "Y", "Z"]
oneQGateindices={"H":0, "P":1, "X":2, "Y":3, "Z":4}


twoQGate_ = ["CNOT", "CZ"]
twoQGateindices={"CNOT":0, "CZ":1}

pauliNoise_ = ["I","X", "Y", "Z"]
pauliNoiseindices={"I":0,"X":1, "Y":2, "Z":3}


class SingeQGate:
    def __init__(self, gateindex, qubitindex):
        self._name = oneQGate_[gateindex]
        self._qubitindex = qubitindex

    def __str__(self):
        return self._name + "[" + str(self._qubitindex) + "]"


class TwoQGate:
    def __init__(self, gateindex, control, target):
        self._name = twoQGate_[gateindex]
        self._control = control
        self._target = target

    def __str__(self):
        return self._name + "[" + str(self._control) + "," + str(self._target)+ "]"


class pauliNoise:
    def __init__(self, noiseindex, qubitindex):
        self._name="n"+str(noiseindex)
        self._noiseindex= noiseindex
        self._qubitindex = qubitindex
        self._noisetype=0


    def set_noisetype(self, noisetype):
        self._noisetype=noisetype


    def __str__(self):
        return self._name +"("+pauliNoise_[self._noisetype] +")" +"[" + str(self._qubitindex) + "]"


class Measurement:
    def __init__(self,measureindex ,qubitindex):
        self._name="M"+str(measureindex)
        self._qubitindex = qubitindex
        self._measureindex=measureindex

    def __str__(self):
        return self._name + "[" + str(self._qubitindex) + "]"


class Reset:
    def __init__(self, qubitindex):
        self._name="R"
        self._qubitindex = qubitindex

    def __str__(self):
        return self._name + "[" + str(self._qubitindex) + "]"



#Class: CliffordCircuit
class CliffordCircuit:


    def __init__(self, qubit_num):
        self._qubit_num = qubit_num
        self._totalnoise=0
        self._totalMeas=0
        self._totalgates=0
        self._gatelists=[]
        self._error_rate=0
        self._index_to_noise={}
        self._index_to_measurement={}

        #self._index_to_measurement={}

        self._shownoise=False
        self._syndromeErrorTable={}
        #Store the repeat match group
        #For example, if we require M0=M1, M2=M3, then the match group is [[0,1],[2,3]]
        self._parityMatchGroup=[]
        self._observable=[]

        self._stim_str=None
        self._stimcircuit=stim.Circuit()


        #self._error_channel


    def set_stim_str(self, stim_str):
        self._stim_str=stim_str


    def set_error_rate(self, error_rate):
        self._error_rate=error_rate

    def get_stim_circuit(self):
        return self._stimcircuit


    def set_observable(self, observablemeasurements):
        self._observable=observablemeasurements


    def get_observable(self):
        return self._observable


    def set_parityMatchGroup(self, parityMatchGroup):
        self._parityMatchGroup=parityMatchGroup

    def get_parityMatchGroup(self):
        return self._parityMatchGroup

    def get_qubit_num(self):
        return self._qubit_num
    
    def get_totalnoise(self):
        return self._totalnoise

    def get_totalMeas(self):
        return self._totalMeas

    '''
    Read the circuit from a file
    Example of the file:

    NumberOfQubit 6
    cnot 1 2
    cnot 1 3
    cnot 1 0
    M 0
    cnot 1 4
    cnot 2 4
    M 4
    cnot 2 5
    cnot 3 5
    M 5
    R 4
    R 5
    cnot 1 4
    cnot 2 4
    M 4
    cnot 2 5
    cnot 3 5
    M 5

    '''
    def read_circuit_from_file(self, filename):
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines
                
                if line.startswith("NumberOfQubit"):
                    # Extract the number of qubits
                    self._qubit_num = int(line.split()[1])
                else:
                    # Parse the gate operation
                    parts = line.split()
                    gate_type = parts[0]
                    qubits = list(map(int, parts[1:]))
                    
                    if gate_type == "cnot":
                        self.add_cnot(qubits[0], qubits[1])
                    elif gate_type == "M":
                        self.add_measurement(qubits[0])
                    elif gate_type == "R":
                        self.add_reset(qubits[0])
                    elif gate_type == "H":
                        self.add_hadamard(qubits[0])
                    elif gate_type == "P":
                        self.add_phase(qubits[0])
                    elif gate_type == "CZ":
                        self.add_cz(qubits[0], qubits[1])
                    elif gate_type == "X":
                        self.add_paulix(qubits[0])
                    elif gate_type == "Y":
                        self.add_pauliy(qubits[0])
                    elif gate_type == "Z":
                        self.add_pauliz(qubits[0])
                    else:
                        raise ValueError(f"Unknown gate type: {gate_type}")

    
    '''
    Compile from a stim circuit string.
    '''
    def compile_from_stim_circuit_str(self, stim_str):
        #self._totalnoise=0
        self._totalnoise=0
        self._totalMeas=0
        self._totalgates=0       

        lines = stim_str.splitlines()
        output_lines = []
        maxum_q_index=0
        '''
        First, read and compute the parity match group and the observable
        '''
        parityMatchGroup=[]
        observable=[]

        
        measure_index_to_line={}
        measure_line_to_measure_index={}             
        current_line_index=0
        current_measure_index=0
        for line in lines:
            stripped_line = line.strip()
            if not stripped_line:
                # Skip empty lines (optional: you could also preserve them)
                current_line_index+=1
                continue
            
            # Keep lines that we do NOT want to split
            if (stripped_line.startswith("TICK") or
                stripped_line.startswith("DETECTOR(") or
                stripped_line.startswith("QUBIT_COORDS(") or                
                stripped_line.startswith("OBSERVABLE_INCLUDE(")):
                current_line_index+=1
                continue

            tokens = stripped_line.split()
            gate = tokens[0]

            if gate == "M":
                measure_index_to_line[current_measure_index]=current_line_index
                measure_line_to_measure_index[current_line_index]=current_measure_index
                current_measure_index+=1

            current_line_index+=1
        

        current_line_index=0
        measure_stack=[]
        for line in lines:
            stripped_line = line.strip()
            if stripped_line.startswith("DETECTOR("):
                meas_index = [token.strip() for token in stripped_line.split() if token.strip().startswith("rec")]
                meas_index = [int(x[4:-1]) for x in meas_index]
                parityMatchGroup.append([measure_line_to_measure_index[measure_stack[x]] for x in meas_index])
                current_line_index+=1
                continue
            elif stripped_line.startswith("OBSERVABLE_INCLUDE("):
                meas_index = [token.strip() for token in stripped_line.split() if token.strip().startswith("rec")]
                meas_index = [int(x[4:-1]) for x in meas_index]
                observable=[measure_line_to_measure_index[measure_stack[x]] for x in meas_index]
                current_line_index+=1
                continue


            tokens = stripped_line.split()
            gate = tokens[0]
            if gate == "M":
                measure_stack.append(current_line_index)
            current_line_index+=1

        '''
        Insert gates
        '''
        for line in lines:
            stripped_line = line.strip()
            if not stripped_line:
                # Skip empty lines (optional: you could also preserve them)
                continue

            # Keep lines that we do NOT want to split
            if (stripped_line.startswith("TICK") or
                stripped_line.startswith("DETECTOR(") or
                stripped_line.startswith("QUBIT_COORDS(") or     
                stripped_line.startswith("OBSERVABLE_INCLUDE(")):
                output_lines.append(stripped_line)
                continue

            tokens = stripped_line.split()
            gate = tokens[0]


            if gate == "CX":
                control = int(tokens[1])
                maxum_q_index=maxum_q_index if maxum_q_index>control else control
                target = int(tokens[2])
                maxum_q_index=maxum_q_index if maxum_q_index>target else target
                self.add_cnot(control, target)


            elif gate == "M":
                qubit = int(tokens[1])
                maxum_q_index=maxum_q_index if maxum_q_index>qubit else qubit
                self.add_measurement(qubit)

            elif gate == "H":
                qubit = int(tokens[1])
                maxum_q_index=maxum_q_index if maxum_q_index>qubit else qubit
                self.add_hadamard(qubit)            

            elif gate == "S":
                qubit = int(tokens[1])
                maxum_q_index=maxum_q_index if maxum_q_index>qubit else qubit
                self.add_phase(qubit)    

            
            elif gate == "R":
                qubits = int(tokens[1])
                maxum_q_index=maxum_q_index if maxum_q_index>qubits else qubits
                self.add_reset(qubits)
            
        '''
        Finally, compiler detector and observable
        '''
        self._parityMatchGroup=parityMatchGroup
        self._observable=observable
        self._qubit_num=maxum_q_index+1
        self.compile_detector_and_observable()    




    def save_circuit_to_file(self, filename):
        pass



    def set_noise_type(self, noiseindex, noisetype):
        self._index_to_noise[noiseindex].set_noisetype(noisetype)


    def reset_noise_type(self):
        for i in range(self._totalnoise):
            self._index_to_noise[i].set_noisetype(0)

    def show_all_noise(self):
        for i in range(self._totalnoise):
            print(self._index_to_noise[i])


    def add_xflip_noise(self, qubit):
        self._stimcircuit.append("X_ERROR", [qubit], self._error_rate)
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1           



    def add_depolarize(self, qubit):
        self._stimcircuit.append("DEPOLARIZE1", [qubit], self._error_rate)
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1        


    def add_cnot_no_noise(self, control, target):
        self._gatelists.append(TwoQGate(twoQGateindices["CNOT"], control, target))
        self._stimcircuit.append("CNOT", [control, target])        



    def add_cnot(self, control, target):
        self._stimcircuit.append("DEPOLARIZE1", [control], self._error_rate)
        self._gatelists.append(pauliNoise(self._totalnoise, control))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1
        self._gatelists.append(pauliNoise(self._totalnoise, target))
        self._stimcircuit.append("DEPOLARIZE1", [target], self._error_rate)
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1
        self._gatelists.append(TwoQGate(twoQGateindices["CNOT"], control, target))
        self._stimcircuit.append("CNOT", [control, target])


    def add_hadamard(self, qubit):
        self._stimcircuit.append("DEPOLARIZE1", [qubit], self._error_rate)        
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1        
        self._gatelists.append(SingeQGate(oneQGateindices["H"], qubit))
        self._stimcircuit.append("H", [qubit])

    def add_phase(self, qubit):
        self._stimcircuit.append("DEPOLARIZE1", [qubit], self._error_rate)   
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1      
        self._gatelists.append(SingeQGate(oneQGateindices["P"], qubit))
        self._stimcircuit.append("S", [qubit])

    def add_cz(self, qubit1, qubit2):
        self._gatelists.append(pauliNoise(self._totalnoise, qubit1))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1
        self._gatelists.append(pauliNoise(self._totalnoise, qubit1))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1
        self._gatelists.append(TwoQGate(twoQGateindices["CZ"], qubit1, qubit2))     


    def add_paulix(self, qubit):
        self._stimcircuit.append("DEPOLARIZE1", [qubit], self._error_rate)   
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1     
        self._gatelists.append(SingeQGate(oneQGateindices["X"], qubit))
        self._stimcircuit.append("X", [qubit])

    def add_pauliy(self, qubit):
        self._stimcircuit.append("DEPOLARIZE1", [qubit], self._error_rate)  
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1    
        self._gatelists.append(SingeQGate(oneQGateindices["Y"], qubit))
        self._stimcircuit.append("Y", [qubit])

    def add_pauliz(self, qubit):
        self._stimcircuit.append("DEPOLARIZE1", [qubit], self._error_rate)  
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1    
        self._gatelists.append(SingeQGate(oneQGateindices["Z"], qubit))
        self._stimcircuit.append("Z", [qubit])


    def add_measurement_no_noise(self, qubit):
        self._gatelists.append(Measurement(self._totalMeas,qubit))
        self._stimcircuit.append("M", [qubit])
        #self._stimcircuit.append("DETECTOR", [stim.target_rec(-1)])
        self._index_to_measurement[self._totalMeas]=self._gatelists[-1]
        self._totalMeas+=1



    def add_measurement(self, qubit):
        self._stimcircuit.append("DEPOLARIZE1", [qubit], self._error_rate)  
        self._gatelists.append(pauliNoise(self._totalnoise, qubit))
        self._index_to_noise[self._totalnoise]=self._gatelists[-1]
        self._totalnoise+=1   
        self._gatelists.append(Measurement(self._totalMeas,qubit))
        self._stimcircuit.append("M", [qubit])
        #self._stimcircuit.append("DETECTOR", [stim.target_rec(-1)])
        self._index_to_measurement[self._totalMeas]=self._gatelists[-1]
        self._totalMeas+=1


    def compile_detector_and_observable(self):
        totalMeas=self._totalMeas
        #print(totalMeas)
        for paritygroup in self._parityMatchGroup:
            #print(paritygroup)
            #print([k-totalMeas for k in paritygroup])
            self._stimcircuit.append("DETECTOR", [stim.target_rec(k-totalMeas) for k in paritygroup])

        self._stimcircuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(k-totalMeas) for k in self._observable], 0)

        #print(self._stimcircuit)


    def add_reset(self, qubit):
        self._gatelists.append(Reset(qubit))
        self._stimcircuit.append("R", [qubit])

    def setShowNoise(self, show):
        self._shownoise=show

    def __str__(self):
        str=""
        for gate in self._gatelists:
            if isinstance(gate, pauliNoise) and not self._shownoise:
                continue
            str+=gate.__str__()+"\n"
        return str


    def get_yquant_latex(self):
        """
        Convert the circuit (stored in self._gatelists) into a yquant LaTeX string.
        This version simply prints each gate (or noise box) in the order they appear,
        without grouping or any fancy logic.
        """
        lines = []
        # Begin the yquant environment
        lines.append("\\begin{yquant}")
        lines.append("")
        
        # Declare qubits and classical bits.
        # Note: Literal braces in the LaTeX code are escaped by doubling them.
        lines.append("% -- Qubits and classical bits --")
        lines.append("qubit {{$\\ket{{q_{{\\idx}}}}$}} q[{}];".format(self._qubit_num))
        lines.append("cbit {{$c_{{\\idx}} = 0$}} c[{}];".format(self._totalMeas))
        lines.append("")
        lines.append("% -- Circuit Operations --")
        
        # Process each gate in the order they were added.
        for gate in self._gatelists:
            if isinstance(gate, pauliNoise):
                # Print the noise box only if noise output is enabled.
                if self._shownoise:
                    lines.append("[fill=red!80]")
                    # The following format string produces, e.g.,:
                    # "box {$n_{8}$} q[2];"
                    lines.append("box {{$n_{{{}}}$}} q[{}];".format(gate._noiseindex, gate._qubitindex))
            elif isinstance(gate, TwoQGate):
                # Two-qubit gate (e.g., CNOT or CZ).
                if gate._name == "CNOT":
                    # Note: yquant syntax for a CNOT is: cnot q[target] | q[control];
                    line = "cnot q[{}] | q[{}];".format(gate._target, gate._control)
                elif gate._name == "CZ":
                    line = "cz q[{}] | q[{}];".format(gate._target, gate._control)
                lines.append(line)
            elif isinstance(gate, SingeQGate):
                # Single-qubit gate.
                if gate._name == "H":
                    line = "h q[{}];".format(gate._qubitindex)

                lines.append(line)
            elif isinstance(gate, Measurement):
                # Measurement is output as three separate lines.
                lines.append("measure q[{}];".format(gate._qubitindex))
                lines.append("cnot c[{}] | q[{}];".format(gate._measureindex, gate._qubitindex))
                lines.append("discard q[{}];".format(gate._qubitindex))
            elif isinstance(gate, Reset):
                # Reset is output as an initialization command.
                lines.append("init {{$\\ket0$}} q[{}];".format(gate._qubitindex))
            else:
                continue
        
        lines.append("")
        lines.append("\\end{yquant}")
        
        return "\n".join(lines)


'''
A
'''
def transpile_stim_with_noise_vector(stimString,noise_vector,totalnoise):
    
    lines = stimString.strip().split('\n')

    MindexToLine={}

    current_line_index=0
    current_M_index=0

    total_detector=0
    for line in lines:
        if line.startswith('M'):
            MindexToLine[current_M_index]=current_line_index
            current_M_index+=1
        if line.startswith('DETECTOR'):     
            total_detector+=1
        current_line_index+=1
 


    s = stim.TableauSimulator(seed=0)
    current_noise_index=0

    detector_result=[]
    observableparity=0
    newstimstr=""
    for line in lines:
        if line.startswith('M'):
            qubit_index=int(line.split(' ')[1])
            s.measure(qubit_index)
            newstimstr+=line+"\n"
        if line.startswith('CX'):
            qubit_index1=int(line.split(' ')[1])
            qubit_index2=int(line.split(' ')[2])
            s.cnot(qubit_index1,qubit_index2)
            newstimstr+=line+"\n"
        if line.startswith('H'):
            qubit_index=int(line.split(' ')[1])
            s.h(qubit_index)
            newstimstr+=line+"\n"
        if line.startswith('S'):
            qubit_index=int(line.split(' ')[1])
            s.s(qubit_index)
            newstimstr+=line+"\n"

        if line.startswith('R'):
            split=line.split(' ')
            for i in range(1,len(split)):
                qubit_index=int(split[i])
                s.reset_z(qubit_index)
            newstimstr+=line+"\n"

        if line.startswith('DETECTOR'):
            split=line.split(' ')        
            parity=0
            for i in range(1,len(split)):
                meas=int(split[i][4:-1])
                tmpmeas=s.current_measurement_record()[meas]
                if tmpmeas:
                    parity+=1
            if parity%2==1:
                detector_result.append(1)
            else:
                detector_result.append(0)

        if line.startswith("OBSERVABLE_INCLUDE(0)"):
            split=line.split(' ')        
            for i in range(1,len(split)):
                meas=int(split[i][4:-1])
                tmpmeas=s.current_measurement_record()[meas]
                if tmpmeas:
                    observableparity+=1
            observableparity=observableparity%2
        if line.startswith('DEPOLARIZE1'):
            split=line.split(' ')
            qubit_index1=int(split[1])
            if noise_vector[current_noise_index]==1:
                s.x(qubit_index1)
                newstimstr+="X error "+str(qubit_index1) +"\n"
            elif noise_vector[current_noise_index+totalnoise]==1:
                s.y(qubit_index1)
                newstimstr+="Y error "+str(qubit_index1) +"\n"
            elif noise_vector[current_noise_index+2*totalnoise]==1:
                s.z(qubit_index1)
                newstimstr+="Z error "+str(qubit_index1) +"\n"
            current_noise_index+=1
            if len(split)==3:
                qubit_index2=int(split[2])
                if noise_vector[current_noise_index]==1:
                    s.x(qubit_index2)
                    newstimstr+="X error "+str(qubit_index2) +"\n"
                elif noise_vector[current_noise_index+totalnoise]==1:
                    s.y(qubit_index2)
                    newstimstr+="Y error "+str(qubit_index2) +"\n"
                elif noise_vector[current_noise_index+2*totalnoise]==1:
                    s.z(qubit_index2)
                    newstimstr+="Z error "+str(qubit_index2) +"\n"
                current_noise_index+=1


    print("-----------------------------New stim circuit:---------------------------------")
    print(newstimstr)

    measurement_result=s.current_measurement_record()


    detector_result.append(observableparity)
    return detector_result




def rewrite_stim_code(code: str) -> str:
    """
    Rewrites a Stim program so that each line contains at most one gate or measurement.
    Lines starting with TICK, R, DETECTOR(, and OBSERVABLE_INCLUDE( are kept as-is.
    Multi-target lines for CX, M, and MR are split up.
    """
    lines = code.splitlines()
    output_lines = []

    for line in lines:
        stripped_line = line.strip()
        if not stripped_line:
            # Skip empty lines (optional: you could also preserve them)
            continue

        # Keep lines that we do NOT want to split
        if (stripped_line.startswith("TICK") or
            stripped_line.startswith("DETECTOR(") or
            stripped_line.startswith("QUBIT_COORDS(") or     
            stripped_line.startswith("OBSERVABLE_INCLUDE(")):
            output_lines.append(stripped_line)
            continue
        
        if (stripped_line.startswith("X_ERROR") or
            stripped_line.startswith("DEPOLARIZE1") or
            stripped_line.startswith("DEPOLARIZE2") or
            stripped_line.startswith("SHIFT_COORDS")            
            ):
            continue
            

        tokens = stripped_line.split()
        gate = tokens[0]

        # Handle 2-qubit gate lines like "CX 0 1 2 3 4 5 ..."
        if gate == "CX":
            qubits = tokens[1:]
            # Pair up the qubits [q0, q1, q2, q3, ...] => (q0,q1), (q2,q3), ...
            for i in range(0, len(qubits), 2):
                q1, q2 = qubits[i], qubits[i + 1]
                output_lines.append(f"CX {q1} {q2}")

        # Handle multi-qubit measurements "M 1 3 5 ..." => each on its own line
        elif gate == "M":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"M {q}")


        elif gate == "MX":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"H {q}")
                output_lines.append(f"M {q}")

        elif gate == "MY":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"S {q}")
                output_lines.append(f"S {q}")
                output_lines.append(f"S {q}")
                output_lines.append(f"H {q}")                
                output_lines.append(f"M {q}")



        elif gate == "H":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"H {q}")

        elif gate == "S":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"S {q}")            

        # Handle multi-qubit measure+reset "MR 1 3 5 ..." => each on its own line
        elif gate == "MR":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"M {q}")
                output_lines.append(f"R {q}")

        elif gate == "R":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"R {q}")
        
        elif gate == "RX":
            qubits = tokens[1:]
            for q in qubits:
                output_lines.append(f"R {q}")
                output_lines.append(f"H {q}")                


        else:
            # If there's some other gate we don't specifically handle,
            # keep it as is, or add more logic if needed.
            output_lines.append(stripped_line)

    return "\n".join(output_lines)



import random

def python_sample_fixed_one_two_three(N, k):
    """
    Returns a list of length N containing exactly k ones 
    (and N-k zeros), in a random order.
    """
    # Step 1: Create a list of k ones and N-k zeros
    arr = [1]*k + [0]*(N-k)
    
    # Step 2: Create a list of 1 or two
    arrtype=[]
    
    for i in range(N):
        arrtype.append(random.randint(1,3))

    
    # Step 2: Shuffle the list randomly
    random.shuffle(arr)
    random.shuffle(arrtype)
    
    return [a * b for a, b in zip(arr, arrtype)]



'''
Test the correctness of detector matrix with Stim simulator

Idea: Sample a noise, convert the circuit to circuit only with these noise, 
and then compare the result of detector with Stim simulator.
'''
def test_with_stim_tableau():
    circuit=CliffordCircuit(4)

    distance=3
    circuit.set_error_rate(0.0001)  
    stim_circuit=stim.Circuit.generated("surface_code:rotated_memory_z",rounds=3*distance,distance=distance).flattened()
    stim_circuit=rewrite_stim_code(str(stim_circuit))
    #print(stim_circuit)
    circuit.set_stim_str(stim_circuit)
    circuit.compile_from_stim_circuit_str(stim_circuit)
    stimcircuit=circuit.get_stim_circuit()
    print("----------------- Original circuit-----------------------------")
    #print(stimcircuit)

    print("Total detectors: ", len(circuit.get_parityMatchGroup())+1)
    print("Total noise: ", 3*circuit.get_totalnoise())

    detectorMatrix=np.array(return_detector_matrix(str(stimcircuit)))
    detectorMatrix = detectorMatrix.T          # or: np.transpose(detector_matrix)
    
    print("Detector matrix: ", detectorMatrix)
    print("Detector matrix shape: ", detectorMatrix.shape)

    '''
    First step, sample a noise
    '''
    totalnoise=circuit.get_totalnoise()
    print("Total noise: ", totalnoise)


    for W in range(1,int(totalnoise/2)):
        for i in range(0,40):
            random_index=python_sample_fixed_one_two_three(totalnoise,W)
            noise_vector=np.array([0]*3*totalnoise)
            for i in range(totalnoise):
                if random_index[i]==1:
                    noise_vector[i]=1
                elif random_index[i]==2:
                    noise_vector[i+totalnoise]=1
                elif random_index[i]==3:
                    noise_vector[i+2*totalnoise]=1    

            print("Noise vector: ", noise_vector)



            detector_result=transpile_stim_with_noise_vector(str(stimcircuit),noise_vector,totalnoise)

            print("-------------Detector result from stim: -------------")
            print(detector_result)


            #print(dectectorMatrix.shape, noise_vector.shape)
            mydetectorresult=np.matmul(detectorMatrix, noise_vector)%2    


            print("-------------My Detector result: -------------")
            print(list(mydetectorresult))

            assert((detector_result==list(mydetectorresult)))




def test_small_circuit(circuit_file_path):
    stim_str=""
    with open(circuit_file_path, "r", encoding="utf-8") as f:
        stim_str = f.read()
    
    circuit=CliffordCircuit(4)     
    circuit.compile_from_stim_circuit_str(stim_str)
    stimcircuit=circuit.get_stim_circuit()
    print("----------------- Original circuit-----------------------------")
    print(stimcircuit)

    print("Total detectors: ", len(circuit.get_parityMatchGroup())+1)
    print("Total noise: ", 3*circuit.get_totalnoise())

    detectorMatrix=np.array(return_detector_matrix(str(stimcircuit)))
    detectorMatrix = detectorMatrix.T          # or: np.transpose(detector_matrix)
    
    print("Detector matrix: ", detectorMatrix)
    print("Detector matrix shape: ", detectorMatrix.shape)

    '''
    First step, sample a noise
    '''
    totalnoise=circuit.get_totalnoise()
    print("Total noise: ", totalnoise)


    for W in range(1,totalnoise):
        for i in range(0,40):
            random_index=python_sample_fixed_one_two_three(totalnoise,W)
            noise_vector=np.array([0]*3*totalnoise)
            for i in range(totalnoise):
                if random_index[i]==1:
                    noise_vector[i]=1
                elif random_index[i]==2:
                    noise_vector[i+totalnoise]=1
                elif random_index[i]==3:
                    noise_vector[i+2*totalnoise]=1    

            print("Noise vector: ", noise_vector)



            detector_result=transpile_stim_with_noise_vector(str(stimcircuit),noise_vector,totalnoise)

            print("-------------Detector result from stim: -------------")
            print(detector_result)


            #print(dectectorMatrix.shape, noise_vector.shape)
            mydetectorresult=np.matmul(detectorMatrix, noise_vector)%2    


            print("-------------My Detector result: -------------")
            print(list(mydetectorresult))

            assert((detector_result==list(mydetectorresult)))    






if __name__ == "__main__":
    #test_small_circuit("C:/Users/yezhu/GitRepos/Sampling/stimprograms/cnoth0")

    test_with_stim_tableau()    