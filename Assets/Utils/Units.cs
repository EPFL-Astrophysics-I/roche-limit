// Class for astronomical units
public static class Units
{
    // Constants
    public static double newtonG_SI = 0.000000000066743; // m^3 / kg / s^2
    public static double au_SI = 149597870700.0;         // m
    public static double earth_to_moon_SI = 384399000.0; // m
    public static double r_sun_SI = 696340000.0;         // m
    public static double m_sun_SI = 1.98847e30;          // kg
    public static double r_earth_SI = 6371000.0;         // m
    public static double m_earth_SI = 5.9722e24;         // kg
    public static double r_moon_SI = 1737400.0;          // m
    public static double m_moon_SI = 7.342e22;           // kg
    public static double year_SI = 31556952.0;           // s
    public static double month_SI = year_SI / 12;        // s
    public static double day_SI = 86400.0;               // s

    // Unit options
    public enum UnitTime { Year, Month, Day }
    public enum UnitLength { AU, SolarRadius, EarthRadius, EarthToMoon }
    public enum UnitMass { SolarMass, EarthMass }

    public static double NewtonG(UnitTime unitTime, UnitLength unitLength, UnitMass unitMass)
    {
        // Time
        double t = year_SI;
        if (unitTime == UnitTime.Month)
        {
            t = month_SI;
        }
        else if (unitTime == UnitTime.Day)
        {
            t = day_SI;
        }

        // Length
        double l = au_SI;
        if (unitLength == UnitLength.SolarRadius)
        {
            l = r_sun_SI;
        }
        else if (unitLength == UnitLength.EarthRadius)
        {
            l = r_earth_SI;
        }
        else if (unitLength == UnitLength.EarthToMoon)
        {
            l = earth_to_moon_SI;
        }

        // Mass
        double m = m_sun_SI;
        if (unitMass == UnitMass.EarthMass)
        {
            m = m_earth_SI;
        }

        return newtonG_SI * m * t * t / l / l / l;
    }

    public static double EarthRotationPeriod(UnitTime unitTime)
    {
        double result = 0;
        switch (unitTime)
        {
            case UnitTime.Day:
                result = 1;
                break;
            case UnitTime.Month:
                result = day_SI / month_SI;
                break;
            case UnitTime.Year:
                result = day_SI / year_SI;
                break;
            default:
                break;
        }
        return result;
    }

    public static double EarthRadius(UnitLength unitLength)
    {
        double result = 0;
        switch (unitLength)
        {
            case UnitLength.AU:
                result = r_earth_SI / au_SI;
                break;
            case UnitLength.EarthRadius:
                result = 1;
                break;
            case UnitLength.SolarRadius:
                result = r_earth_SI / r_sun_SI;
                break;
            default:
                break;
        }
        return result;
    }

    public static double EarthMass(UnitMass unitMass)
    {
        double result = 1;
        if (unitMass == UnitMass.SolarMass)
        {
            result = m_earth_SI / m_sun_SI;
        }
        return result;
    }

    public static double LunarRadius(UnitLength unitLength)
    {
        double result = 0;
        switch (unitLength)
        {
            case UnitLength.AU:
                result = r_moon_SI / au_SI;
                break;
            case UnitLength.EarthRadius:
                result = r_moon_SI / r_earth_SI;
                break;
            case UnitLength.SolarRadius:
                result = r_moon_SI / r_sun_SI;
                break;
            default:
                break;
        }
        return result;
    }

    public static double LunarMass(UnitMass unitMass)
    {
        double result = 0;
        switch (unitMass)
        {
            case UnitMass.EarthMass:
                result = m_moon_SI / m_earth_SI;
                break;
            case UnitMass.SolarMass:
                result = m_moon_SI / m_sun_SI;
                break;
            default:
                break;
        }
        return result;
    }

    public static double LunarDistance(UnitLength unitLength)
    {
        double result = 0;
        switch (unitLength)
        {
            case UnitLength.AU:
                result = earth_to_moon_SI / au_SI;
                break;
            case UnitLength.EarthRadius:
                result = earth_to_moon_SI / r_earth_SI;
                break;
            case UnitLength.SolarRadius:
                result = earth_to_moon_SI / r_sun_SI;
                break;
            default:
                break;
        }
        return result;
    }
}
