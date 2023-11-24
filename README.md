# Quite
⚛️ Quick and Dirty Quantum Computer Emulator

Heavily based on https://github.com/davidbkemp/jsqubits

## Example Usage

Grover Algorithm searching for `x² - 8x + 12 = 0`
```py
from quite import Q

inputbits = range(1,9)

def amplify(qstate, f):
	return (qstate
		.apply(inputbits, 0, f)
		.h(inputbits)
		.apply(inputbits, 0, lambda x: x==0)
		.h(inputbits))

def search_f(x):
	return x**2 - 8*x + 12 == 0

q = Q('000000001').h(True)

for _ in range(10):
	for i in range(100):
		measure, _ = q[inputbits]
		if search_f(measure):
			print(f'found {measure} in {i} attempts')
			break
		q = amplify(q, search_f)
```
