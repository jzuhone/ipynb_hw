import numpy as np
from JSAnimation import IPython_display
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import animation
font = {'family' : 'serif',
        'size'   : 14}
matplotlib.rc('font', **font)

# Constants

sec_per_day = 24*3600.
sec_per_year = 365.25*sec_per_day
mass_sun_cgs = 1.989e33
G = 6.6726e-8
cm_per_au = 1.49597871e13
onethird = 1./3.
onesixth = 1./6.

def ellipse_area(x, e, p):
    if e == 1.:
        return 0.5*p*np.tan(0.5*x)*(onesixth*np.cos(0.5*x)**(-2)+onethird)
    f = np.abs(1.-e)*np.tan(0.5*x)/np.sqrt(np.abs(1.-e*e))
    if e < 1.:
        f = np.arctan(f)
    else:
        f = -np.arctanh(f)
    f /= (np.abs(1.-e*e))**1.5
    g = 0.5*e*np.sin(x)/((e*e-1.)*(e*np.cos(x)+1.))
    return p*(f+g)

class KeplersLaws(object):

    def __init__(self):
        pass

    def problem1(self, e=None):
        self.e = e
        self.draw_areas = False
        self._compute_orbit()
        return self._draw_plot()
    
    def problem2(self, e=None, t1=None, t2=None, t3=None, t4=None):
        self.e = e
        self.t1 = t1
        self.t2 = t2
        self.t3 = t3
        self.t4 = t4
        self.draw_areas = True
        self._compute_orbit()
        return self._draw_plot()

    def problem3(self, e=None):
        self.e = e
        self.draw_areas = False
        self._compute_orbit()
        print "a = %g AU, P = %g years" % (self.a/cm_per_au, self.tmax/sec_per_year)
        
    def _compute_orbit(self):
        x = cm_per_au
        y = 0.0
        vx = 0.0
        mu = G*mass_sun_cgs
        vy = np.sqrt(mu*(self.e+1.)/x)
        self.p = (x*vy)**2/mu/cm_per_au
        if self.e < 1.:
            self.a = x/(1.-self.e)
            ra = self.a*(1.+self.e)
            self.tmax = 2.*np.pi*self.a*np.sqrt(self.a/mu)
            self.lim = 1.5*ra/cm_per_au
        else:
            self.tmax = 10*sec_per_year
            self.lim = 10*x/cm_per_au
        t = 0.0
        X = [x]
        Y = [y]
        VX = [vx]
        VY = [vy]
        T = [t]

        first = True

        gx_old = 0.0
        gy_old = 0.0

        while t <= self.tmax:
            r = np.sqrt(x*x+y*y)
            g = -mu/(r*r)
            gx = g*x/r
            gy = g*y/r
            dt = 0.001*np.sqrt(-r/g)
            if first:
                wterm = 0.5*dt
                woldterm = 0.0
                first = False
            else:
                wterm = 0.5*dt + onethird*dt_old + onesixth*dt**2/dt_old
                woldterm = onesixth*(dt_old**2-dt**2)/dt_old
            vx += wterm*gx + woldterm*gx_old
            vy += wterm*gy + woldterm*gy_old
            x += vx*dt
            y += vy*dt
            gx_old = gx
            gy_old = gy
            dt_old = dt
            X.append(x)
            Y.append(y)
            VX.append(vx)
            VY.append(vy)
            t += dt
            T.append(t)
        
        self.T = np.array(T)/sec_per_year
        self.X = np.array(X)/cm_per_au
        self.Y = np.array(Y)/cm_per_au
        self.R = np.sqrt(self.X*self.X+self.Y*self.Y)

        self.tt = np.linspace(0, self.tmax, 100)/sec_per_year
        self.xx = np.interp(self.tt, self.T, self.X)
        self.yy = np.interp(self.tt, self.T, self.Y)
        self.rr = np.sqrt(self.xx*self.xx+self.yy*self.yy)

    def _draw_plot(self):
        fig = plt.figure(figsize=(8,6))
        self.ax = plt.axes(xlim=(-self.lim, self.lim), ylim=(-self.lim, self.lim))
        self.sun = plt.Circle((0.0,0.0), radius=0.25, facecolor='Yellow', edgecolor='k')
        self.ax.add_patch(self.sun)
        self.line, = self.ax.plot([], [], lw=2, color='Green')
        self.ax.set_aspect('equal')
        self.ax.set_xlabel('x (AU)')
        self.ax.set_ylabel('y (AU)')
        self.radius_txt = self.ax.text(0.55, 0.1, "", transform=self.ax.transAxes)
        self.time_txt = self.ax.text(0.55, 0.2, "", transform=self.ax.transAxes)
        self.earth = plt.Circle((0.0,0.0), radius=0.1, facecolor='Blue', edgecolor='k')
        if self.draw_areas:
            self.line1 = plt.Line2D((0,self.earth.center[0]),(0,self.earth.center[1]), lw=2, color="r")
            self.line2 = plt.Line2D((0,self.earth.center[0]),(0,self.earth.center[1]), lw=2, color="r")
            self.line3 = plt.Line2D((0,self.earth.center[0]),(0,self.earth.center[1]), lw=2, color="b")
            self.line4 = plt.Line2D((0,self.earth.center[0]),(0,self.earth.center[1]), lw=2, color="b")
            self.coord1 = np.array([np.interp(self.t1, self.T, self.X), np.interp(self.t1, self.T, self.Y)])
            self.coord2 = np.array([np.interp(self.t2, self.T, self.X), np.interp(self.t2, self.T, self.Y)])
            self.coord3 = np.array([np.interp(self.t3, self.T, self.X), np.interp(self.t3, self.T, self.Y)])
            self.coord4 = np.array([np.interp(self.t4, self.T, self.X), np.interp(self.t4, self.T, self.Y)])
            theta1 = np.arctan2(self.coord1[1], self.coord1[0]) 
            theta2 = np.arctan2(self.coord2[1], self.coord2[0]) 
            theta3 = np.arctan2(self.coord3[1], self.coord3[0]) 
            theta4 = np.arctan2(self.coord4[1], self.coord4[0]) 
            if theta1 > 0.0 and theta2 < 0.0:
                self.area1 = 2*ellipse_area(np.pi, self.e, self.p)-ellipse_area(theta2, self.e, self.p)-ellipse_area(theta1, self.e, self.p)
            else:
                self.area1 = ellipse_area(theta2, self.e, self.p)-ellipse_area(theta1, self.e, self.p)
            if theta3 > 0.0 and theta4 < 0.0:
                self.area2 = 2*ellipse_area(np.pi, self.e, self.p)+ellipse_area(theta4, self.e, self.p)-ellipse_area(theta3, self.e, self.p)
            else:
                self.area2 = ellipse_area(theta4, self.e, self.p)-ellipse_area(theta3, self.e, self.p)

        return animation.FuncAnimation(fig, self._animate, init_func=self._init,
                                       frames=len(self.tt), interval=50) 

    def _init(self):
        self.line.set_data([], [])
        self.radius_txt.set_text("")
        self.time_txt.set_text("")
        self.earth.center = (1.0, 0.0)
        self.ax.add_artist(self.earth)
        return (self.line, self.time_txt, self.radius_txt)
    
    def _animate(self, i):
        self.line.set_data(self.xx[:i+1], self.yy[:i+1])
        self.earth.center = (self.xx[i], self.yy[i])
        self.ax.add_artist(self.earth)
        if self.draw_areas:
            if self.tt[i] >= self.t1:
                self.line1.set_xdata((0.0,self.coord1[0]))
                self.line1.set_ydata((0.0,self.coord1[1]))
                self.ax.add_artist(self.line1)
            if self.tt[i] >= self.t1 and self.tt[i] <= self.t2:
                self.line2.set_xdata((0.0,self.earth.center[0]))
                self.line2.set_ydata((0.0,self.earth.center[1]))
                self.ax.add_artist(self.line2)
            if self.tt[i] > self.t2:
                a1 = self.ax.text(0.1, 0.2, r"$\mathrm{A_1\ =\ %g}$" % (self.area1), color="r", transform=self.ax.transAxes)
                self.ax.add_artist(a1)
            if self.tt[i] >= self.t3:
                self.line3.set_xdata((0.0,self.coord3[0]))
                self.line3.set_ydata((0.0,self.coord3[1]))
                self.ax.add_artist(self.line3)
            if self.tt[i] >= self.t3 and self.tt[i] <= self.t4:
                self.line4.set_xdata((0.0,self.earth.center[0]))
                self.line4.set_ydata((0.0,self.earth.center[1]))
                self.ax.add_artist(self.line4)
            if self.tt[i] > self.t4:
                a2 = self.ax.text(0.1, 0.1, r"$\mathrm{A_2\ =\ %g}$" % (self.area2), color="b", transform=self.ax.transAxes)
                self.ax.add_artist(a2)
        self.radius_txt.set_text("r = %5.2f AU" % (self.rr[i]))
        self.time_txt.set_text("t = %5.2f years" % (self.tt[i]))
        return (self.line, self.time_txt, self.radius_txt)


