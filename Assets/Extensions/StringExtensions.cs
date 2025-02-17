using System.Linq;
using System.Text;

using System;

/// <summary>
/// A collection of extension methods to obfuscate a string. These methods are applied
/// (in a fixed order) to the activeLettersFilter string so that you can type a sentence
/// and have it obfuscated before mapping its letters to sphere points.
/// </summary>
public static class StringExtensions
{
    /// <summary>
    /// Removes vowels from the string.
    /// </summary>
    public static string Disemvowel(this string str)
    {
        return new string(str.Where(c => !"aeiouAEIOU".Contains(c)).ToArray());
    }

    /// <summary>
    /// Removes consecutive duplicate characters.
    /// </summary>
    public static string Squeeze(this string str)
    {
        StringBuilder sb = new StringBuilder();
        foreach (char c in str)
        {
            if (sb.Length == 0 || sb[sb.Length - 1] != c)
                sb.Append(c);
        }
        return sb.ToString();
    }

    /// <summary>
    /// Removes duplicate characters (only the first occurrence remains).
    /// </summary>
    public static string Deduplicate(this string str)
    {
        return new string(str.Distinct().ToArray());
    }

    /// <summary>
    /// Removes all whitespace from the string.
    /// </summary>
    public static string RemoveWhitespace(this string str)
    {
        StringBuilder sb = new StringBuilder();
        foreach (char c in str)
        {
            if (!char.IsWhiteSpace(c))
                sb.Append(c);
        }
        return sb.ToString();
    }

    /// <summary>
    /// Randomly shuffles the characters in the string.
    /// </summary>
    public static string Rearrange(this string str)
    {
        char[] chars = str.ToCharArray();
        System.Random rng = new System.Random();
        for (int i = chars.Length - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (chars[i], chars[j]) = (chars[j], chars[i]);
        }
        return new string(chars);
    }

    /// <summary>
    /// Shifts each letter by a given offset (or random if not provided).
    /// </summary>
    public static string Shift(this string str, int? offset = null)
    {
        int actualOffset = offset ?? new System.Random().Next(26);
        StringBuilder sb = new StringBuilder();
        foreach (char c in str)
        {
            if (char.IsLetter(c))
            {
                char baseChar = char.IsUpper(c) ? 'A' : 'a';
                int shifted = (c - baseChar + actualOffset) % 26;
                sb.Append((char)(baseChar + shifted));
            }
            else
                sb.Append(c);
        }
        return sb.ToString();
    }

    /// <summary>
    /// Reverses the string.
    /// </summary>
    public static string Reverse(this string str)
    {
        char[] arr = str.ToCharArray();
        Array.Reverse(arr);
        return new string(arr);
    }

    /// <summary>
    /// Randomly shuffles the order of words in the string.
    /// </summary>
    public static string CutUp(this string str)
    {
        string[] words = str.Split(' ');
        System.Random rng = new System.Random();
        for (int i = words.Length - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (words[i], words[j]) = (words[j], words[i]);
        }
        return string.Join(" ", words);
    }

    /// <summary>
    /// Rotates the words in the string by a given offset (or random if not provided).
    /// </summary>
    public static string Rotate(this string str, int? offset = null)
    {
        string[] words = str.Split(' ');
        int actualOffset = offset ?? new System.Random().Next(words.Length);
        string[] rotated = new string[words.Length];
        for (int i = 0; i < words.Length; i++)
            rotated[(i + actualOffset) % words.Length] = words[i];
        return string.Join(" ", rotated);
    }

    /// <summary>
    /// Reverses the order of words in the string.
    /// </summary>
    public static string Flip(this string str)
    {
        string[] words = str.Split(' ');
        Array.Reverse(words);
        return string.Join(" ", words);
    }
}