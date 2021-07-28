                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            if unit == 'kB':
                val = int(val, 10)
                val /= 1024.
                unit = 'MB'
                val = '%.0f' % val
            txt.append('%s: %s %s' % (k, val, unit))
        return ', '.join([] + txt)

class IoMeas(object):
    '''
    A measurement class for the *Time* class that measures disk I/O.
    '''    
    def __init__(self):
        self.io0 = get_procio()
    def format_diff(self, other):
        txt = []
        d1 = self.io0
        d0 = other.io0
        for k,knice in [('rchar',None), ('wchar',None),
                        ('read_bytes','rb'), ('write_bytes', 'wb')]:
            v1 = d1.get(k)
            v0 = d0.get(k)
            unit = 'b'
            dv = float(v1 - v0)
            for uu in ['kB', 'MB', 'GB']:
                if dv < 2048:
                    break
                dv = dv / 1024.
                unit = uu
            val = '%.3g' % dv
            if knice is None:
                kk = k
            else:
                kk = knice
            txt.append('%s: %s %s' % (kk, val, unit))
        return ', '.join([] + txt)

class CpuMeas(object):
    '''
    A measurement class for the *Time* class that measures CPU and Wall tim.
    '''    
    def __init__(self):
        import datetime
        from time import process_time
        self.wall = datetime.datetime.now()
        self.cpu = process_time()

    def cpu_seconds_since(self, other):
        '''
        Returns the CPU time used since the *other* timestamp was taken.
        eg::
            t0 = CpuMeas()
            # do stuff
            t1 = CpuMeas()
            dt = t1.cpu_seconds_since(t0)
            print('That took', dt, 'seconds of CPU time')
        '''
        return self.cpu - other.cpu

    def wall_seconds_since(self, other):
        '''
        Returns the wall-clock time in seconds since the *other* timestamp
        was taken.
        eg::
        
            t0 = CpuMeas()
            # do stuff
            t1 = CpuMeas()
            dt = t1.wall_seconds_since(t0)
            print('That took', dt, 'seconds')
        '''
        dwall = (self.wall - other.wall)
        # python2.7
        if hasattr(dwall, 'total_seconds'):
            dwall = dwall.total_seconds()
        else:
            dwall = (dwall.microseconds + (dwall.seconds + dwall.days * 24. * 3600.) * 1e6) / 1e6
        return dwall
        
    def format_diff(self, other):
        dwall = self.wall_seconds_since(other)
        dcpu = self.cpu_seconds_since(other)
        return 'Wall: %.2f s, CPU: %.2f s' % (dwall, dcpu)
        
class Time(object):
    '''
    A class for recording how much time (or other resources) are used
    by a process.  The *CpuMeas* class is used by default; others can
    be added as desired.
    Use like this::
        from astrometry.util.ttime import Time, MemMeas
        # measure memory usage too
        Time.add_measurement(MemMeas)
        
        def mymethod():
            t0 = Time()
            # do stuff
            t1 = Time()
            print('Time taken:', t1-t0)
    '''
    @staticmethod
    def add_measurement(m):
        Time.measurements.append(m)
    @staticmethod
    def remove_measurement(m):
        Time.measurements.remove(m)
    measurements = [CpuMeas]

    def __init__(self):
        self.meas = [m() for m in Time.measurements]

    def __sub__(self, other):
        '''
        Returns a string representation of the difference: self - other
        '''
        meas = ', '.join([m.format_diff(om) for m,om in zip(self.meas, other.meas)])
        return meas
© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
