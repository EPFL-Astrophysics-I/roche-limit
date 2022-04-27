using TMPro;
using UnityEngine;
using UnityEngine.EventSystems;

public class LanguageToggle : MonoBehaviour, IPointerClickHandler
{
    [SerializeField] private TextMeshProUGUI labelEN = default;
    [SerializeField] private TextMeshProUGUI labelFR = default;

    [Header("Light Theme")]
    [SerializeField] private Color lightActiveColor = Color.black;
    [SerializeField] private Color lightInactiveColor = Color.gray;

    [Header("Dark Theme")]
    [SerializeField] private Color darkActiveColor = Color.gray;
    [SerializeField] private Color darkInactiveColor = Color.black;

    public enum ActiveLanguage { EN, FR }
    private ActiveLanguage activeLanguage = ActiveLanguage.EN;
    private bool canToggle = true;

    private Color active;
    private Color inactive;

    public void OnPointerClick(PointerEventData eventData)
    {
        if (canToggle)
        {
            ToggleLanguageDisplay();
            SendMessageUpwards("SetLanguage", activeLanguage.ToString(), SendMessageOptions.DontRequireReceiver);
        }
    }

    private void Awake()
    {
        if (!labelEN)
        {
            Debug.LogWarning("Need to assign EN TMP");
            canToggle = false;
        }

        if (!labelFR)
        {
            Debug.LogWarning("Need to assign FR TMP");
            canToggle = false;
        }
    }

    private void ToggleLanguageDisplay()
    {
        if (activeLanguage == ActiveLanguage.EN)
        {
            labelEN.color = inactive;
            labelFR.color = active;
            activeLanguage = ActiveLanguage.FR;
        }
        else
        {
            labelEN.color = active;
            labelFR.color = inactive;
            activeLanguage = ActiveLanguage.EN;
        }
    }

    private void SetLanguageDisplay()
    {
        if (activeLanguage == ActiveLanguage.EN)
        {
            labelEN.color = active;
            labelFR.color = inactive;
        }
        else
        {
            labelEN.color = inactive;
            labelFR.color = active;
        }
    }

    public void SetLightTheme()
    {
        active = lightActiveColor;
        inactive = lightInactiveColor;
        SetLanguageDisplay();
    }

    public void SetDarkTheme()
    {
        active = darkActiveColor;
        inactive = darkInactiveColor;
        SetLanguageDisplay();
    }
}
