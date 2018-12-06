package rnaSeqGUI;

public final class FlagForCodeSubstitution {
	/*
	 * Creates an object that describes a flag and a value to replace the flag. A
	 * flag is a String that will be replaced by the String valueToSubstitute. This
	 * is performed by the applyAllFlagSubstitutions method, which uses the
	 * String.replace method to replace the flag with the valueToSubstitute
	 */
	private final String flag;
	private final String valueToSubstitute;

	public FlagForCodeSubstitution(String flag, String valueToSubstitute) {
		this.flag = flag;
		this.valueToSubstitute = valueToSubstitute;
	}

	public String getFlag() {
		return flag;
	}

	public String getValueToSubstitute() {
		return valueToSubstitute;
	}

}
