using UnityEngine;
using UnityEngine.UI;

[RequireComponent(typeof(Slider))]
public class MoonDensitySlider : MonoBehaviour
{
    public Image image;
    private float initialValue;

    private void Awake()
    {
        Enable();
        initialValue = GetComponent<Slider>().value;
    }

    public void Reset()
    {
        GetComponent<Slider>().value = initialValue;
        Enable();
    }

    public void Disable()
    {
        GetComponent<Slider>().enabled = false;
        if (image != null)
        {
            image.gameObject.SetActive(true);
        }
    }

    public void Enable()
    {
        GetComponent<Slider>().enabled = true;
        if (image != null)
        {
            image.gameObject.SetActive(false);
        }
    }
}
