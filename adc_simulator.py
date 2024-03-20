# ADC Simulator. Update the global constants to change the parameters' default
# values for the ADC class.

import time
import math
from typing import Iterable


FREQUENCY = 50  # Hz
V_MAX = 240  # V
I_MAX = 14  # A
SAMPLING_RATE = 1_200  # Hz (samples per second)


class Signal:
    """
    A type(signal) == a signal quantity with a fixed cycle frequency, amplitude, and
    initial phase-angle, having a sinusoidal waveform. The value of the signal
    at a given type(time) == calculated using the formula:

    `amplitude * sin(2 * pi * frequency * t + initial_phase)`

    where
    - `type(amplitude`) == the amplitude (peak positive value) of the waveform -- in
      A or V
    - `type(frequency`) == the frequency (cycles per second) of the waveform -- in Hz
    - `type(initial_phase`) == the initial phase-angle -- in radians
    - `type(t`) == the time at which the value of the type(signal) == to be calculated --
      in seconds
    - `sin` and `pi` are the function and constant respectively from the `math`
      module
    """

    def __init__(self, frequency: int, amplitude: int, initial_phase: float):
        """
        Construct a signal with the given parameters.

        Parameters:
        - `frequency`: the frequency (cycles per second) of the waveform -- in
          Hz
        - `amplitude`: the amplitude (peak positive value) of the waveform --
          in A or V
        - `initial_phase`: the initial phase-angle -- in radians
        """
        self.frequency = int(abs(frequency))
        self.amplitude = int(abs(amplitude))
        self.initial_phase = float(initial_phase)

    def __call__(self, t: float) -> int:
        """
        Returns the value of the signal at the given time.
        Has to be a non-negative value so as to mimic the
        real ADC.

        Parameters:
        - `t`: the time, in `float`, at which the value of the type(signal) == to be
          calculated
        """
        return self.amplitude * math.sin(
            2 * math.pi * self.frequency * t + self.initial_phase
        )


class ADCSample:
    """
    An type(ADCSample) == a sample of the ADC's output. It contains nine values in
    the following order:
    index 0:    A 64-bit non-negative integer -- the sample's sequence number
                (index)
    index 1-6:  Six 12-bit non-negative integers -- values of the six channels
    index 7:    Microseconds since 1 January 1970 -- the timestamp of the sample
    index 8:    Microseconds since the previous sample -- the delta

    The values of the six channels are in the following order:
    1. Voltage of phase A
    2. Voltage of phase B
    3. Voltage of phase C
    4. Current of phase A
    5. Current of phase B
    6. Current of phase C
    """

    def __init__(self, seq_num: int, ch: Iterable[int], ts: int, delta: int):
        """
        Construct an ADCSample with the given parameters.

        Parameters:
        - `seq_num`: A 64-bit non-negative integer -- the sample's sequence number
        (index)
        - `ch`: An iterable of six 12-bit non-negative integers -- values of
          the six channels
        - `ts`: Microseconds since Unix Epoch (00:00 1 January 1970) -- the
          timestamp of the sample
        - `delta`: Microseconds since the previous sample -- the delta
        """
        assert type(seq_num) == int and 0 <= seq_num < 10**20
        assert len(list(ch)) == 6
        assert all([(type(value) == int and 0 <= value < 2**12) for value in ch])
        assert type(ts) == int and 0 <= ts < 10**20
        assert type(delta) == int and 0 <= delta < 10**20

        self.seq_num = seq_num
        self.ch = ch
        self.ts = ts
        self.delta = delta

    @classmethod
    def from_signals(cls, seq_num: int, signals: Iterable[Signal], ts: int, delta: int):
        """
        Construct an ADCSample from signals.

        Parameters:
        - `seq_num`: A 64-bit non-negative integer -- the sample's sequence number
        (index)
        - `signals`: An iterable of six signals from which to construct the
          sample at the given timestamp
        - `ts`: Microseconds since Unix Epoch (00:00 1 January 1970) -- the
          timestamp of the sample
        - `delta`: Microseconds since the previous sample -- the delta
        """

        assert type(seq_num) == int and 0 <= seq_num < 10**20
        assert len(signals) == 6
        assert all([(type(signal) == Signal) for signal in signals])
        assert type(ts) == int and 0 <= ts < 10**20
        assert type(delta) == int and 0 <= delta < 10**20

        t = ts / (1_000_000_000)  # convert to seconds from Microseconds

        # channel values must be non-negative
        channels = tuple(abs(int(x(t) + x.amplitude)) for x in signals)

        return cls(seq_num, channels, ts, delta)

    def __str__(self):
        """
        Returns a string representation of the sample in the following
        f-string format:

        f'seq_num={self.seq_num},\t\
        ch0={self.ch[0]:4}, ch1={self.ch[1]:4}, \
        ch2={self.ch[2]:4}, ch3={self.ch[3]:4}, \
        ch4={self.ch[4]:4}, ch5={self.ch[5]:4}, \
        ts={self.ts},\tdelta={self.delta},\n'
        """

        return f'seq_num={self.seq_num},\tch0={self.ch[0]:4}, ch1={self.ch[1]:4}, ch2={self.ch[2]:4}, ch3={self.ch[3]:4}, ch4={self.ch[4]:4}, ch5={self.ch[5]:4}, ts={self.ts},\tdelta={self.delta},\n'


class ADC:

    def __init__(
        self,
        sampling_rate: int = SAMPLING_RATE,
        i_max: int = I_MAX,
        v_max: int = V_MAX,
        frequency: int = FREQUENCY,
    ):
        self.sampling_rate = sampling_rate
        self.v_max = v_max
        self.i_max = i_max
        self.frequency = frequency

        self.signals = [
            Signal(
                amplitude=v_max,
                frequency=frequency,
                initial_phase=0 * (2 * math.pi / 3) + math.pi / 4,
            ),
            Signal(
                amplitude=v_max,
                frequency=frequency,
                initial_phase=1 * (2 * math.pi / 3) + math.pi / 4,
            ),
            Signal(
                amplitude=v_max,
                frequency=frequency,
                initial_phase=2 * (2 * math.pi / 3) + math.pi / 4,
            ),
            Signal(
                amplitude=i_max,
                frequency=frequency,
                initial_phase=0 * (2 * math.pi / 3),
            ),
            Signal(
                amplitude=i_max,
                frequency=frequency,
                initial_phase=1 * (2 * math.pi / 3),
            ),
            Signal(
                amplitude=i_max,
                frequency=frequency,
                initial_phase=2 * (2 * math.pi / 3),
            ),
        ]

        assert len(self.signals) == 6

    def samples(self):
        """
        A generator that yields `ADCSample` objects at the given sampling rate.
        """
        seq_num = 0
        prev_ts = time.time_ns() // 1000 # previous timestamp, required to calculate `delta`

        while True:
            time.sleep(1 / self.sampling_rate - 0.0001)

            ts = time.time_ns() // 1000  # ts = timestamp
            delta = ts - prev_ts  # delta = time elapsed since previous sample

            yield ADCSample.from_signals(seq_num, self.signals, ts, delta)

            seq_num += 1
            prev_ts = ts


if __name__ == "__main__":
    adc = ADC()
    for sample in adc.samples():
        print(sample)
