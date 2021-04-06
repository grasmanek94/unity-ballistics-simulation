using UnityEngine;
using System;

public class BallisticTrajectory : MonoBehaviour
{
    // http://mackila.com/airsoft/atp/01-d-01.htm

    public GameObject looker;
    public float mass = 0.00020f; // kg
    public float diameter = 0.006f; // m
    public float initialSpeed = 100.584f; // m/s, 330 FPS = 100.584 m/s
    public float avgDragCoefficient = 0.47f;
    public float pitch = 45.0f; // Z
    private float gravity = -9.80655f; // m/s^2 TODO: Replace with GetGravity
    private float area; //m^2
    public float currentSpeedX = 0.0f;
    public float currentSpeedY = 0.0f;
    public float time;

    float GetAirDensity(float altitude)
    {
        float L = -0.0065f; // K/m
        float T0 = 288.15f; // K
        float M = 0.028964f; // kg/mol
        float R = 8.31447f; // J/(mol * K)
        float p0 = 101325.0f; // kg/(ms^2)

        float T = T0 + L * altitude; // K

        float p = p0 * Mathf.Pow(1.0f + ((L * altitude) / T0), (-gravity * M) / (R * -L));

        float rho = (p * M) / (R * T);

        return rho; // kg/m^3
    }

    float GetForceOfDrag(float air_density, float speed)
    {
        return 0.5f * avgDragCoefficient * air_density * area * speed * speed;
    }

    float GetAirViscosity(float temperature)
    {
        float C1 = 1.458e-6f;
        float C2 = 110.4f;
        return (C1 * Mathf.Pow(temperature, 1.5f)) / (temperature + C2);
    }

    float ReynoldsNumber(float diameter, float air_density, float speed, float air_viscosity)
    {
        return diameter * air_density * speed * air_viscosity;
    }

    float GetVelocity(float initial_speed, float acceleration, float time)
    {
        return initial_speed + acceleration * time;
    }

    float GetDistance(float initial_distance, float velocity_a, float time, float acceleration)
    {
        return initial_distance + velocity_a * time + 0.5f * acceleration * time * time;
    }

    private float DegreeToRadian(float angle)
    {
        return Mathf.PI * angle / 180.0f;
    }

    float GetGravity(float latitude)
    {
        float rlat = DegreeToRadian(latitude);

        float sl = 0.005278895f * Mathf.Pow(Mathf.Sin(rlat), 2.0f);
        float ssl = 0.0000589f * Mathf.Pow(Mathf.Sin(2.0f * rlat), 2.0f);

        return 9.7803185f * (1.0f + (sl - ssl));
    }

    float GetMagnusForce(float air_density, float speed)
    {
        float CL = -0.0020907f;
        return CL * air_density * speed * speed * area;
    }

    float GetTerminalVelocity(float air_density)
    {
        double part_a = mass * gravity;
        double part_b = 0.5f * avgDragCoefficient * air_density * area;

        return -(float)Math.Sqrt(Math.Abs(part_a / part_b));
    }

    // Start is called before the first frame update
    void Start()
    {
        float rad_pitch = DegreeToRadian(pitch);

        currentSpeedX = initialSpeed * Mathf.Cos(rad_pitch);
        currentSpeedY = initialSpeed * Mathf.Sin(rad_pitch);

        area = Mathf.Pow(diameter / 2.0f, 2.0f) * Mathf.PI;
        time = 0.0f;

        looker.transform.LookAt(transform);
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        if(transform.position.y < 0.0f)
        {
            return;
        }

        time += Time.fixedDeltaTime;

        float delta = Time.deltaTime;
        float altitude = 300.0f;

        float air_density = GetAirDensity(altitude);

        float terminalVelocity = GetTerminalVelocity(air_density);

        if(currentSpeedY < terminalVelocity)
        {
            currentSpeedY = terminalVelocity;
        }

        float speed = Mathf.Sqrt(currentSpeedX * currentSpeedX + currentSpeedY * currentSpeedY);
        
        float dragX = GetForceOfDrag(air_density, currentSpeedX);
        float dragY = GetForceOfDrag(air_density, currentSpeedY);

        // float magnus = GetMagnusForce(air_density, speed);
        // float magnusX = magnus * Mathf.Abs(currentSpeedY) / speed;
        // float magnusY = magnus * Mathf.Abs(currentSpeedX) / speed;

        float gravityForce = mass * gravity;

        float totalXForce = 
            // magnusX +
            -dragX;

        float totalYForce = 
            gravityForce +
            // -magnusY +
            -dragY;

        float accelerationX = totalXForce / mass;
        float accelerationY = totalYForce / mass;

        currentSpeedX += accelerationX * delta;
        currentSpeedY += accelerationY * delta;

        if (currentSpeedY < terminalVelocity)
        {
            currentSpeedY = terminalVelocity;
        }

        transform.localPosition = new Vector3(
            transform.localPosition.x + delta * currentSpeedX,
            transform.localPosition.y + delta * currentSpeedY,
            0.0f
       );

       looker.transform.LookAt(transform);
    }
}
