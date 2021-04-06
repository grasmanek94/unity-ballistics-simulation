using UnityEngine;
using System;

public class BallisticTrajectory3D : MonoBehaviour
{
    // http://mackila.com/airsoft/atp/01-d-01.htm

    public GameObject looker;
    public double mass = 0.00020; // kg
    public double diameter = 0.006; // m
    public double initialSpeed = 100.584; // m/s, 330 FPS = 100.584 m/s
    public double avgDragCoefficient = 0.47;
    public double yaw = 0.0; // direction in degrees
    public double pitch = 45.0; // pitch in degrees
    public double windSpeed = 0.0; // m/s
    public double windDirection = 0.0; // degrees (which direction the wind flows to relative from yaw)
    public double temperature = 18.0; // degrees C
    public double altitude = 300.0; // meters above sea level
    private double gravity = -9.80655; // m/s^2 TODO: Replace with GetGravity
    private double area; //m^2
    private double windSpeedX;
    private double windSpeedZ;
    public double currentSpeedX = 0.0;
    public double currentSpeedY = 0.0;
    public double currentSpeedZ = 0.0;
    public double time;

    double GetAirDensity(double altitude, double temperature)
    {
        double L = -0.0065; // K/m
        double T0 = 273.15 + temperature; // K
        double M = 0.028964; // kg/mol
        double R = 8.31447; // J/(mol * K)
        double p0 = 101325.0; // kg/(ms^2)

        double T = T0 + L * altitude; // K

        double p = p0 * Math.Pow(1.0 + ((L * altitude) / T0), (-gravity * M) / (R * -L));

        double rho = (p * M) / (R * T);

        return rho; // kg/m^3
    }

    double GetForceOfDrag(double air_density, double speed)
    {
        return Math.Sign(speed) * 0.5 * avgDragCoefficient * air_density * area * speed * speed;
    }

    double GetAirViscosity(double temperature)
    {
        double C1 = 1.458e-6;
        double C2 = 110.4;
        return (C1 * Math.Pow(temperature, 1.5)) / (temperature + C2);
    }

    double ReynoldsNumber(double diameter, double air_density, double speed, double air_viscosity)
    {
        return diameter * air_density * speed * air_viscosity;
    }

    double GetVelocity(double initial_speed, double acceleration, double time)
    {
        return initial_speed + acceleration * time;
    }

    double GetDistance(double initial_distance, double velocity_a, double time, double acceleration)
    {
        return initial_distance + velocity_a * time + 0.5 * acceleration * time * time;
    }

    private double DegreeToRadian(double angle)
    {
        return Mathf.PI * angle / 180.0;
    }

    double GetGravity(double latitude)
    {
        double rlat = DegreeToRadian(latitude);

        double sl = 0.005278895f * Math.Pow(Math.Sin(rlat), 2.0);
        double ssl = 0.0000589f * Math.Pow(Math.Sin(2.0 * rlat), 2.0);

        return 9.7803185f * (1.0 + (sl - ssl));
    }

    double GetMagnusForce(double air_density, double speed)
    {
        double CL = -0.0020907;
        return CL * air_density * speed * speed * area;
    }

    double GetTerminalVelocity(double air_density)
    {
        double part_a = mass * gravity;
        double part_b = 0.5 * avgDragCoefficient * air_density * area;

        return -Math.Sqrt(Math.Abs(part_a / part_b));
    }

    // Start is called before the first frame update
    void Start()
    {
        double rad_pitch = DegreeToRadian(pitch);
        double rad_yaw = DegreeToRadian(yaw);

        double currentForwardSpeed = initialSpeed * Math.Cos(rad_pitch);
        currentSpeedX = currentForwardSpeed * Math.Cos(rad_yaw);
        currentSpeedY = initialSpeed * Math.Sin(rad_pitch);
        currentSpeedZ = currentForwardSpeed * Math.Sin(rad_yaw);

        double rad_wind = DegreeToRadian(windDirection);

        windSpeedX = windSpeed * Math.Cos(rad_wind);
        windSpeedZ = windSpeed * Math.Sin(rad_wind);

        area = Math.Pow(diameter / 2.0, 2.0) * Math.PI;
        time = 0.0;

        looker.transform.LookAt(transform);
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        if (transform.position.y < 0.0f)
        {
            return;
        }

        time += Time.fixedDeltaTime;

        double delta = Time.deltaTime;

        double air_density = GetAirDensity(altitude, temperature);
        double terminalVelocity = GetTerminalVelocity(air_density);

        if (currentSpeedY < terminalVelocity)
        {
            currentSpeedY = terminalVelocity;
        }

        double speed = Math.Sqrt(
            currentSpeedX * currentSpeedX + 
            currentSpeedY * currentSpeedY +
            currentSpeedZ * currentSpeedZ
        );
 
        double dragX = -GetForceOfDrag(air_density, currentSpeedX);
        double dragY = -GetForceOfDrag(air_density, currentSpeedY);
        double dragZ = -GetForceOfDrag(air_density, currentSpeedZ);

        double windDragX = GetForceOfDrag(air_density, windSpeedX);
        double windDragZ = GetForceOfDrag(air_density, windSpeedZ);

        double speedXZ = Math.Sqrt(
            currentSpeedX * currentSpeedX +
            currentSpeedZ * currentSpeedZ
        );


        // this is wrong.. but oh well
        double magnus = GetMagnusForce(air_density, speed);

        double magnusXZ = magnus * Math.Abs(currentSpeedY) / speed;
        double magnusY = magnus * Math.Abs(speedXZ) / speed;
        double magnusX = magnusXZ * currentSpeedX / speedXZ;
        double magnusZ = magnusXZ * currentSpeedZ / speedXZ;
        // what we actually need to to above is, to know the spin rate and axis, 
        // and apply the force perpendicular to the spin axis
        // those forces need to be recalculated/realigned to the XYZ speed vectors,
        // that's a lot of work that I don't want to do atm
        // (not to mention I don't have a way to measure the spin).

        double gravityForce = mass * gravity;

        double totalXForce =
            magnusX +
            windDragX +
            dragX;

        double totalYForce =
            gravityForce +
            -magnusY +
            dragY;

        double totalZForce =
            magnusZ +
            windDragZ +
            dragZ;

        double accelerationX = totalXForce / mass;
        double accelerationY = totalYForce / mass;
        double accelerationZ = totalZForce / mass;

        currentSpeedX += accelerationX * delta;
        currentSpeedY += accelerationY * delta;
        currentSpeedZ += accelerationZ * delta;

        if (currentSpeedY < terminalVelocity)
        {
            currentSpeedY = terminalVelocity;
        }

        transform.localPosition = new Vector3(
            transform.localPosition.x + (float)(delta * currentSpeedX),
            transform.localPosition.y + (float)(delta * currentSpeedY),
            transform.localPosition.z + (float)(delta * currentSpeedZ)
       );

        looker.transform.LookAt(transform);
    }
}
