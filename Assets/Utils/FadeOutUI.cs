using System.Collections;
using UnityEngine;

[RequireComponent(typeof(CanvasGroup))]
public class FadeOutUI : MonoBehaviour
{
    public float fadeInTime = 1;
    public float fadeInDelay = 0;
    public float fadeOutTime = 1;
    public float fadeOutDelay = 0;

    private CanvasGroup canvasGroup;
    private bool doneFading;

    private void OnEnable()
    {
        if (!canvasGroup)
        {
            canvasGroup = GetComponent<CanvasGroup>();
        }

        Reset();
    }

    public void TriggerFadeOut()
    {
        StopAllCoroutines();
        StartCoroutine(FadeOut(fadeOutTime, fadeOutDelay));
    }

    public void TriggerFadeIn()
    {
        if (canvasGroup.alpha != 1)
        {
            StartCoroutine(FadeIn(fadeInTime, fadeInDelay));
        }
    }

    public void TriggerReset(float delay = 0)
    {
        if (canvasGroup.alpha != 1)
        {
            Invoke(nameof(Reset), delay);
        }
    }

    public void SetAlpha(float alpha)
    {
        StopAllCoroutines();
        canvasGroup.alpha = alpha;
    }

    private void Reset()
    {
        StopAllCoroutines();
        canvasGroup.alpha = 0;
    }

    private IEnumerator FadeOut(float fadeTime, float startDelay)
    {
        yield return new WaitForSeconds(startDelay);

        float time = 0;
        float startAlpha = canvasGroup.alpha;

        while (time < fadeTime)
        {
            time += Time.deltaTime;
            canvasGroup.alpha = Mathf.Lerp(startAlpha, 0, time / fadeTime);
            yield return null;
        }

        canvasGroup.alpha = 0;
    }

    private IEnumerator FadeIn(float fadeTime, float startDelay)
    {
        yield return new WaitForSeconds(startDelay);

        float time = 0;
        float startAlpha = canvasGroup.alpha;

        while (time < fadeTime)
        {
            time += Time.deltaTime;
            canvasGroup.alpha = Mathf.Lerp(startAlpha, 1, time / fadeTime);
            yield return null;
        }

        canvasGroup.alpha = 1;
    }
}
