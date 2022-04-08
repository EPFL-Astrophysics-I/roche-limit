using UnityEngine;
using UnityEngine.UI;

public class MoonRingSlideController : SimulationSlideController
{
    private MoonRing sim;
    private MoonRingPrefabs prefabs;

    [SerializeField] private bool useLights = true;
    [SerializeField] private Slider moonDensitySlider;

    private void Awake()
    {
        sim = (MoonRing)simulation;
        if (!sim.TryGetComponent(out prefabs))
        {
            Debug.LogWarning("Could not find MoonRingPrefabs component");
        }
    }

    private void Start()
    {
        if (prefabs && moonDensitySlider)
        {
            float rocheLimit = 1.26f * Mathf.Pow(5.5f / moonDensitySlider.value, 0.33f) * prefabs.earthRadius;
            prefabs.DrawRocheLimit(rocheLimit);
        }
    }

    public override void InitializeSlide()
    {
        if (!prefabs)
        {
            return;
        }

        if (prefabs && moonDensitySlider)
        {
            float rocheLimit = 1.26f * Mathf.Pow(5.5f / moonDensitySlider.value, 0.33f) * prefabs.earthRadius;
            prefabs.DrawRocheLimit(rocheLimit);
        }

        prefabs.SetLightsVisibility(useLights);
    }

    public void SetMoonDensity(float value)
    {
        float rocheLimit = 1.26f * Mathf.Pow(5.5f / value, 0.33f) * prefabs.earthRadius;
        prefabs.DrawRocheLimit(rocheLimit);
        if (rocheLimit > sim.lunarDistance)
        {
            if (moonDensitySlider)
            {
                moonDensitySlider.GetComponent<MoonDensitySlider>().Disable();
            }

            sim.BreakUpMoon();
        }
    }
}
