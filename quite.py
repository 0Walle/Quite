from math import sqrt, cos, sin, log2
import random

def Q(val):
	if type(val) is str:
		return QState.from_string(val)
	if val == 0:
		return QState(1, [complex(1), complex(0)])

	size = int(log2(val))+1
	amps = _make_amps(size)
	amps[val] = complex(1)
	return QState(size, amps)

def _format_complex(c):
	if c.imag:
		return f"({c.real:.4}+{c.imag:.4}i)"
	return f"{c.real:.4g}"

def _make_amps(size):
	return [complex(0)]*(2**size)

def _control_bitmask(bits):
	if not bits: return None
	mask = 0
	for b in bits:
		mask += 1 << b
	return mask


SQRT1_2 = complex(1/sqrt(2))
HADAMARD = [
	[SQRT1_2, SQRT1_2],
	[SQRT1_2, -SQRT1_2]
]

XMATRIX = [
	[complex(0), complex(1)],
	[complex(1), complex(0)]
]

YMATRIX = [
	[complex(0), -1j],
	[1j, complex(0)]
]

ZMATRIX = [
	[complex(1), complex(0)],
	[complex(0), complex(-1)]
]

SMATRIX = [
	[complex(1), complex(0)],
	[complex(0), 1j]
]

TMATRIX = [
	[complex(1), complex(0)],
	[complex(0), complex(SQRT1_2, SQRT1_2)]
]

class QState:
	def from_string(bitstr):
		val = int(bitstr, 2)
		size = len(bitstr)
		amp = [complex(0,0)]*(2**size)
		amp[val] = complex(1,0)
		return QState(size, amp)

	def apply_one(qstate, target, op, control=[]):
		amps = _make_amps(qstate.size)
		bitmask = 1 << target
		cmask = _control_bitmask(control)
		for i, amp in qstate._each():
			if cmask is None or cmask & i == cmask:
				bit = int(i & bitmask > 0)
				zero = i & ~bitmask
				one = i | bitmask
				amps[zero] += op[0][bit] * amp
				amps[one] += op[1][bit] * amp
			else:
				amps[i] = amp
		return QState(qstate.size, amps)

	def __init__(self, size, amp):
		self.size = size
		self.amp = amp

	def _each(self):
		return enumerate(self.amp)

	def norm(self):
		scale = 1 / sqrt(sum(abs(amp) ** 2 for amp in self.amp))
		amps = [ scale*amp for amp in self.amp ]
		return QState(self.size, amps)

	def __str__(self):
		return ' + '.join( f"{_format_complex(amp) if amp != 1 else ''}|{i:0{self.size}b}>" for i, amp in enumerate(self.amp) if amp != 0j)

	__repr__=__str__

	def __add__(self, other):
		amp = [a+b for a, b in zip(self.amp, other.amp)]
		return QState(self.size, amp).norm()

	def __sub__(self, other):
		return self + other * -1

	def __mul__(self, ammount):
		amp = [amp*ammount for amp in self.amp]
		return QState(self.size, amp)

	def __neg__(self):
		return self * -1

	__rmul__ = __mul__

	def __matmul__(self, other):
		amps = _make_amps(self.size+other.size)
		for i, ampa in self.norm()._each():
			for offset, j, ampb in other.norm()._each():
				amps[(i << other.size) + j] = ampa * ampb
		return QState(self.size+other.size, amps)

	def apply_matrix(self, target, op, control=[]):
		if target is True:
			target = range(self.size)

		if type(control) is int: control = [control]
		if type(target) is int: target = [target]

		result = self
		for bit in target:
			result = QState.apply_one(result, bit, op, control)
		return result

	def h(self, t, c=[]):
		return self.apply_matrix(t, HADAMARD, c)

	def x(self, t, c=[]):
		return self.apply_matrix(t, XMATRIX, c)

	def y(self, t, c=[]):
		return self.apply_matrix(t, YMATRIX, c)

	def z(self, t, c=[]):
		return self.apply_matrix(t, ZMATRIX, c)

	def s(self, t, c=[]):
		return self.apply_matrix(t, SMATRIX, c)

	def t(self, t, c=[]):
		return self.apply_matrix(t, TMATRIX, c)

	def r(self, t, angle, c=[]):
		RMATRIX = [
			[complex(1), complex(0)],
			[complex(0), complex(cos(angle), sin(angle))]
		]
		return self.apply_matrix(t, RMATRIX, c)

	def apply(self, input, output, fn):
		if type(input) is int:
			input = range(input, input)

		input_mask = (1 << (input.stop + 1)) - 1
		skip = set()

		amps = _make_amps(self.size)
		for i, amp in self._each():
			if i in skip: continue
			val = (i & input_mask) >> input.start
			result = i ^ ((fn(val) << output) & (1 << output))
			if i == result:
				amps[i] = amp
			else:
				skip.add(result)
				amps[i] = self.amp[result]
				amps[result] = amp

		return QState(self.size, amps)

	def random_state(self):
		rand = random.random()
		acc = 0
		for i, amp in self._each():
			acc += abs(amp)**2
			if acc > rand: break
		return i

	def measure(self, bits=[]):
		rand = self.random_state()
		mask = _control_bitmask(bits)
		mask_rand = rand & mask

		amps = _make_amps(self.size)
		for i, amp in self._each():
			if i & mask == mask_rand:
				amps[i] = amp

		value = 0
		for i in range(self.size-1, -1, -1):
			if i in bits:
				value <<= 1
				value += int(rand & (1 << i) > 0)


		return value, QState(self.size, amps).norm()

	def __getitem__(self, k):
		if type(k) is int:
			return self.measure([k])

		if type(k) is slice:
			return self.measure(range(k.start, k.stop))

		return self.measure(k)

	__int__=random_state

	def __and__(self, k):
		if type(k) is int:
			return self.z(k)
		return self.z(k[0], k[1:])

	def __xor__(self, k):
		return self.x(k)