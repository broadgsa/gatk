package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

public enum EventType {
    BASE_SUBSTITUTION(0, "M"),
    BASE_INSERTION(1, "I"),
    BASE_DELETION(2, "D");

    public final int index;
    private final String representation;

    private EventType(int index, String representation) {
        this.index = index;
        this.representation = representation;
    }

    public static EventType eventFrom(int index) {
        switch (index) {
            case 0:
                return BASE_SUBSTITUTION;
            case 1:
                return BASE_INSERTION;
            case 2:
                return BASE_DELETION;
            default:
                throw new ReviewedStingException(String.format("Event %d does not exist.", index));
        }        
    }
    
    public static EventType eventFrom(String event) {
        for (EventType eventType : EventType.values())
            if (eventType.representation.equals(event))
                return eventType;

        throw new ReviewedStingException(String.format("Event %s does not exist.", event));
    }

    @Override
    public String toString() {
        return representation;
    }
}